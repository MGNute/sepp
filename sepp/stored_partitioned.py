

from sepp.algorithm import AbstractAlgorithm
from sepp.config import options
from sepp.tree import PhylogeneticTree
from sepp.alignment import MutableAlignment, ExtendedAlignment
from sepp.problem import SeppProblem, Problem
from sepp.filemgr import get_temp_file
import dendropy
from dendropy.dataobject.tree import Tree
from sepp.jobs import HMMBuildJob, HMMSearchJob, HMMAlignJob, PplacerJob,\
    MergeJsonJob
from sepp.scheduler import JobPool, Join
from sepp import get_logger
from sepp.math_utils import lcm
from sepp.exhaustive import JoinSearchJobs, JoinAlignJobs, ExhaustiveAlgorithm
import pickle, os

_LOG = get_logger(__name__)

# TODO: make this an option
sepp_inventory_path = '/projects/tallis/nute/work/metagenomics/silva_1231/class_align/sepp_root_inventory.txt'
storage_path = '/projects/tallis/nute/work/metagenomics/silva_1231/class_align/stored_partitioned_sepp'

class PartitionedSeppAlgorithm(ExhaustiveAlgorithm):

    def __init__(self):
        ExhaustiveAlgorithm.__init__(self)
        self.root_job_list={}
        self.parse_inventory()

    def check_options(self, supply = []):
        self.check_outputprefix()


    def parse_inventory(self):
        f=open(sepp_inventory_path,'r')
        jobs = f.readlines()
        f.close()
        for i in jobs:
            a=i.strip().split('\t')
            args = {
                'name':a[0],
                'alignment': a[1],
                'tree': a[2],
                'info': a[3]
            }
            self.root_job_list[a[0]]=args.copy()

    def read_and_divide_fragments(self, chunks, extra_frags = {}):
        _LOG.debug("start reading fragment files and breaking to chunks: %d" %chunks)
        self.root_problem.fragments = MutableAlignment()
        self.root_problem.fragments.read_file_object(self.options.fragment_file)
        for (k,v) in extra_frags.iteritems():
            self.root_problem.fragments[k] = v.replace("-","")
        alg_chunks = self.root_problem.fragments.divide_to_equal_chunks(chunks)
        ret = []
        for i in xrange(0,chunks):
            temp_file = None
            if alg_chunks[i]:
                temp_file = get_temp_file("fragment_chunk_%d" %i, "fragment_chunks", ".fasta")
                alg_chunks[i].write_to_path(temp_file)
            ret.append(temp_file)
        _LOG.debug("fragment files read and divided.")
        return ret

    def read_alignment_and_tree(self, alnfile=None, treefile=None):
        _LOG.info("Reading input alignment: %s" % (alnfile))
        alignment = MutableAlignment()
        alignment.read_file_object(alnfile)

        # fragments = MutableAlignment()
        # fragments.read_file_object(self.options.fragment_file);
        _LOG.info("Reading input tree: %s" % (treefile))
        tree = PhylogeneticTree(dendropy.Tree(stream=open(treefile,'r'),
                                              schema="newick",
                                              preserve_underscores=True))

        return (alignment, tree)

    def _create_root_problem(self, tree, alignment, name):
        ''' Create the root problem'''
        # root_problem = SeppProblem(tree.leaf_node_names(),parent=self.root_problem)
        root_problem = SeppProblem(tree.leaf_node_names())
        root_problem.label = name + '_root'
        root_problem.subalignment = alignment
        root_problem.subtree = tree
        root_problem.major_root_reference=self.root_problem
        return root_problem

    def build_subproblems(self):
        placement_tree_map={}
        self.root_problem = MultiPartSeppRootProblem()
        for part in self.root_job_list.values():
            (alignment, tree) = self.read_alignment_and_tree(part['alignment'], part['tree'])
            assert isinstance(tree, PhylogeneticTree)
            assert isinstance(alignment, MutableAlignment)

            tree.get_tree().resolve_polytomies()
            # Label edges with numbers so that we could assemble things back
            # at the end
            tree.lable_edges()

            ''' Make sure size values are set, and are meaningful. '''
            self.check_and_set_sizes(alignment.get_num_taxa())

            rp = self._create_root_problem(tree, alignment, part['name'])
            self.root_problem.add_child(rp)
            part['problem']=rp

            _LOG.info("Decomposing tree: %s" % (part['name']))
            ''' Decompose the tree based on placement subsets'''
            ptm = PhylogeneticTree(Tree(tree.den_tree), name=part['name']).decompose_tree(
                                            self.options.placement_size,
                                            strategy=self.strategy,
                                            minSize = self.minsubsetsize,
                                            tree_map = {},
                                            named = True)
            # print placement_tree_map.keys()
            assert len(ptm) > 0, ("Tree could not be decomposed"
                    " given the following settings; strategy:%s minsubsetsize:%s placement_size:%s"
                    %(self.strategy, self.minsubsetsize, self.options.placement_size))
            _LOG.info("Breaking into %d placement subsets." %len(ptm))
            placement_tree_map.update(ptm)
            # self.root_problem.add_list_of_children(self.root_job_list)

            ''' For placement subsets create a placement subproblem, and decompose further'''
            for (p_key,p_tree) in ptm.items():
                assert isinstance(p_tree, PhylogeneticTree)
                # placement_problem  = SeppProblem(p_tree.leaf_node_names(), self.root_problem)
                placement_problem  = SeppProblem(p_tree.leaf_node_names(), part['problem'])
                placement_problem.subtree = p_tree
                placement_problem.label = "P_%s" %str(p_key)
                _LOG.debug("Placement subset %s has %d nodes" %(placement_problem.label,len(p_tree.leaf_node_names())))
                ''' Further decompose to alignment subsets '''
                alignment_tree_map = PhylogeneticTree(Tree(p_tree.den_tree)).decompose_tree(
                                            self.options.alignment_size,
                                            strategy=self.strategy,
                                            minSize = self.minsubsetsize,
                                            tree_map = {}, decomp_strategy = self.options.decomp_strategy)
                assert len(alignment_tree_map) > 0, ("Tree could not be decomposed"
                " given the following settings; strategy:%s minsubsetsize:%s alignmet_size:%s"
                %(self.strategy, self.minsubsetsize, self.options.alignment_size))

                _LOG.debug("Placement subset %s has %d alignment subsets: %s" %(placement_problem.label,len(alignment_tree_map.keys()),str(sorted(alignment_tree_map.keys()))))
                _LOG.debug("Placement subset %s has %d taxa:" %(placement_problem.label,sum([len(a_tree.leaf_node_names()) for a_tree in alignment_tree_map.values()])))
                for (a_key, a_tree) in alignment_tree_map.items():
                    assert isinstance(a_tree, PhylogeneticTree)
                    self.modify_tree(a_tree)
                    alignment_problem  = SeppProblem(a_tree.leaf_node_names(),
                                                      placement_problem)
                    alignment_problem.subtree = a_tree
                    alignment_problem.label = "A_%s_%s" %(str(p_key),str(a_key))

        ''' Divide fragments into chunks, to help achieve better parallelism'''
        fragment_chunk_files = self.create_fragment_files()
        for alignment_problem in self.root_problem.iter_leaves():
            for afc in xrange(0,len(fragment_chunk_files)):
                frag_chunk_problem  = SeppProblem(alignment_problem.taxa,
                                              alignment_problem)
                frag_chunk_problem.subtree = alignment_problem.subtree
                frag_chunk_problem.label = alignment_problem.label.replace("A_", "FC_") + "_" +str(afc)
                frag_chunk_problem.fragments = fragment_chunk_files[afc]

        _LOG.info("Breaking into %d alignment subsets." %(len(list(self.root_problem.iter_leaves()))))
        _LOG.info("Breaking each alignment subset into %d fragment chunks." %len(fragment_chunk_files))
        # _LOG.info("Subproblem structure: %s" %str(self.root_problem))
        return self.root_problem

    def build_jobs(self):
        for part in self.root_job_list.values():
            assert isinstance(part['problem'], SeppProblem)
            for placement_problem in part['problem'].get_children():
                ''' Create pplacer jobs'''
                pj = PplacerJob()
                placement_problem.add_job(pj.job_type,pj)
                pj.partial_setup_for_subproblem(placement_problem, part['info'])

                '''For each alignment subproblem, ...'''
                for alg_problem in placement_problem.children:
                    assert isinstance(alg_problem, SeppProblem)
                    ''' create the build model job'''
                    bj = HMMBuildJob()
                    bj.setup_for_subproblem(alg_problem,symfrac=self.symfrac, molecule=self.molecule,**vars(self.options.hmmbuild))
                    alg_problem.add_job(bj.job_type, bj)
                    ''' create the search jobs'''
                    for fc_problem in alg_problem.get_children():
                        sj = HMMSearchJob()
                        sj.partial_setup_for_subproblem(fc_problem.fragments, fc_problem, self.elim, self.filters)
                        fc_problem.add_job(sj.job_type, sj)
                        ''' create the align job'''
                        aj = HMMAlignJob()
                        fc_problem.add_job(aj.job_type, aj)
                        aj.partial_setup_for_subproblem(fc_problem, molecule=self.molecule)

    def _get_new_Join_Align_Job(self):
        return MultiPartJoinAlignJobs()

    def connect_jobs(self):
        ''' a callback function called after hmmbuild jobs are finished'''
        _LOG.info('Connecting jobs...')
        def enq_job_searchfragment(result, search_job):
            search_job.hmmmodel = result
            JobPool().enqueue_job(search_job)
        for part in self.root_job_list.values():
            assert isinstance(part['problem'], SeppProblem)
            for placement_problem in part['problem'].get_children():
                '''For each alignment subproblem, ...'''
                for alg_problem in placement_problem.children:
                    assert isinstance(alg_problem, SeppProblem)
                    ''' create the build model job'''
                    bj = alg_problem.jobs["hmmbuild"]
                    ''' create the search jobs'''
                    for fc_problem in alg_problem.get_children():
                        sj = fc_problem.jobs["hmmsearch"]
                        ''' connect bulid and search jobs'''
                        bj.add_call_Back(lambda result, next_job = sj: enq_job_searchfragment(result, next_job))
                '''Join all align jobs of a placement subset (enqueues placement job)'''
                jaj = self._get_new_Join_Align_Job()
                _LOG.debug('placement problem: %s' % placement_problem.label) #debug
                jaj.setup_with_placement_problem(placement_problem)
            ''' Join all search jobs together (enqueues align jobs)'''
        jsj = MultiPartJoinSearchJobs()
        jsj.setup_with_root_problem(self.root_problem)

    def enqueue_firstlevel_job(self):
        for rp in self.root_problem.children:
            for p in rp.children:
                for ap in p.children:
                    JobPool().enqueue_job(ap.jobs["hmmbuild"])

    def merge_results(self):
        '''

        '''
        _LOG.info('Merging results...')
        self.results = {}
        # print self.root_problem.root_job_list #DEBUG
        for part in self.root_job_list.values():
            rp = part['problem']
            id = part['name']
            _LOG.info('merging ID: %s' % id)
            assert isinstance(rp,SeppProblem)

            '''Generate single extended alignment'''
            fullExtendedAlignment = rp.get_children()[0].jobs["pplacer"].get_attribute("full_extended_alignment_object")
            for pp in rp.get_children()[1:]:
                extended_alignment = pp.jobs["pplacer"].get_attribute("full_extended_alignment_object")
                fullExtendedAlignment.merge_in(extended_alignment,convert_to_string=True)
            self.results[id]=fullExtendedAlignment

            mergeinput = []
            '''Append main tree to merge input'''
            mergeinput.append("%s;" %(rp.subtree.compose_newick(labels = True)))
            jsons = []
            for pp in rp.get_children():
                assert isinstance(pp,SeppProblem)
                if (pp.get_job_result_by_name("pplacer") is None):
                  continue
                '''Append subset trees and json locations to merge input'''
                mergeinput.append("%s;\n%s" %(pp.subtree.compose_newick(labels = True),
                                  pp.get_job_result_by_name("pplacer")))
            mergeinput.append("")
            mergeinput.append("")
            meregeinputstring = "\n".join(mergeinput)
            mergeJsonJob = MergeJsonJob()
            mergeJsonJob.setup(meregeinputstring,
                               self.get_output_filename(id + "_placement.json"))
            mergeJsonJob.run()

    def output_results(self):
        ''' Merged json file is already saved in merge_results function and
            full extended alignment already created in merge_results function
        '''
        _LOG.info('outputing results:')
        for i in self.results.keys():
            _LOG.info(i)
            outfilename = self.get_output_filename(i + "_alignment.fasta")
            self.results[i].write_to_path(outfilename)
            self.results[i].remove_insertion_columns()
            outfilename = self.get_output_filename(i + "_alignment_masked.fasta")
            self.results[i].write_to_path(outfilename)


class MultiPartSeppRootProblem(Problem):
    root_job_list = {}
    def __init__(self, taxa=None, parent=None, root_job_list=None):
        Problem.__init__(self,parent)
        self.label = 'major_root'
        if root_job_list is not None:
            self.add_list_of_children(root_job_list)

    def add_list_of_children(self, root_job_list):
        for i in root_job_list.values():
            self.add_child(i['problem'])

class MultiPartJoinAlignJobs(JoinAlignJobs):
    def __init__(self):
        JoinAlignJobs.__init__(self)

    def merge_subalignments(self):
        '''
        Merge alignment subset extended alignments to get one extended alignment
        for current placement subset.
        '''
        pp = self.placement_problem
        _LOG.info("Merging sub-alignments for placement problem : %s." %(pp.label))
        ''' First assign fragments to the placement problem'''
        pp.fragments = pp.parent.major_root_reference.fragments.get_soft_sub_alignment([])
        frags = []
        for ap in pp.get_children():
            frags.extend(ap.fragments)
        pp.fragments.seq_names.update(frags)
        ''' Then Build an extended alignment by merging all hmmalign results'''
        extendedAlignment = ExtendedAlignment(pp.fragments.seq_names)
        for ap in pp.children:
            assert isinstance(ap, SeppProblem)
            ''' Get all fragment chunk alignments for this alignment subset'''
            aligned_files = [fp.get_job_result_by_name('hmmalign') for
                                fp in ap.children if
                                fp.get_job_result_by_name('hmmalign') is not None]
            _LOG.info("Merging fragment chunks for subalignment : %s." %(ap.label))
            ap_alg = ap.read_extendend_alignment_and_relabel_columns\
                        (ap.jobs["hmmbuild"].infile , aligned_files)
            _LOG.info("Merging alignment subset into placement subset: %s." %(ap.label))
            extendedAlignment.merge_in(ap_alg,convert_to_string=False)
            del ap_alg

        extendedAlignment.from_bytearray_to_string()
        return extendedAlignment

    def setup_with_placement_problem(self, placement_problem):
        self.placement_problem = placement_problem
        for p in placement_problem.iter_leaves():
            # try:
            self.add_job(p.jobs["hmmalign"])
            # except:
            #     print p.jobs.keys()
            #     import sys
            #     sys.exit(1)

class MultiPartJoinSearchJobs(JoinSearchJobs):
    def __init__(self):
        JoinSearchJobs.__init__(self)

    def setup_with_root_problem(self, root_problem):
        self.root_problem = root_problem
        for rp in self.root_problem.get_children():
            for p in rp.iter_leaves():
                self.add_job(p.jobs["hmmsearch"])



    def perform(self):
        '''
        Distributes fragments to alignments subsets with best score,
        and runs align jobs on those. Also, creates new chunks of fragments
        for better parallelism.
        '''

        ''' Figure out which fragment should go to which subproblem'''
        self.figureout_fragment_subset()

        ''' For each alignment subproblem,
        1) make sure its fragments are evenly distributed to fragment chunks.
        2) Setup alignment jobs for its children and enqueue them'''
        alg_problems = [alg for rp in self.root_problem.children for p in rp.children for alg in p.children ]
        for alg_problem in alg_problems:
            assert isinstance(alg_problem, SeppProblem)
            chunks = len(alg_problem.get_children())
            fragment_chunks = alg_problem.fragments.divide_to_equal_chunks(chunks)

            ''' Now setup alignment jobs and enqueue them'''
            for (i,fragment_chunk_problem) in enumerate(alg_problem.children):
                fragment_chunk_problem.fragments = fragment_chunks[i]
                aj = fragment_chunk_problem.jobs['hmmalign']
                assert isinstance(aj,HMMAlignJob)
                ''' First Complete setting up alignments'''
                aj.hmmmodel = alg_problem.get_job_result_by_name('hmmbuild')
                aj.base_alignment = alg_problem.jobs["hmmbuild"].infile

                if fragment_chunk_problem.fragments is None or fragment_chunk_problem.fragments.is_empty():
                    aj.fake_run = True
                else:
                    fragment_chunk_problem.fragments.write_to_path(aj.fragments)
                ''' Now the align job can be put on the queue '''
                JobPool().enqueue_job(aj)

def augment_parser():
    from sepp.config import get_parser
    parser = get_parser()
    parser.description = "This script runs SEPP over a series of large alignment-tree combinations, and allows the target alignment and tree to be decomposed and stored in advance. ."






class StoredPartitionedSeppAlgorithm(PartitionedSeppAlgorithm):
    def __init__(self):

        PartitionedSeppAlgorithm.__init__(self)
        self.open_hmm_job_log()
        self.parse_inventory()

    def build_subproblems_from_scratch(self):
        placement_tree_map={}
        self.root_problem = MultiPartSeppRootProblem()
        for part in self.root_job_list.values():
            (alignment, tree) = self.read_alignment_and_tree(part['alignment'], part['tree'])
            assert isinstance(tree, PhylogeneticTree)
            assert isinstance(alignment, MutableAlignment)

            tree.get_tree().resolve_polytomies()
            # Label edges with numbers so that we could assemble things back
            # at the end
            tree.lable_edges()

            ''' Make sure size values are set, and are meaningful. '''
            self.check_and_set_sizes(alignment.get_num_taxa())

            rp = self._create_root_problem(tree, alignment, part['name'])
            self.root_problem.add_child(rp)
            part['problem']=rp

            _LOG.info("Decomposing tree: %s" % (part['name']))
            ''' Decompose the tree based on placement subsets'''
            ptm = PhylogeneticTree(Tree(tree.den_tree), name=part['name']).decompose_tree(
                                            self.options.placement_size,
                                            strategy=self.strategy,
                                            minSize = self.minsubsetsize,
                                            tree_map = {},
                                            named = True)
            # print placement_tree_map.keys()
            assert len(ptm) > 0, ("Tree could not be decomposed"
                    " given the following settings; strategy:%s minsubsetsize:%s placement_size:%s"
                    %(self.strategy, self.minsubsetsize, self.options.placement_size))
            _LOG.info("Breaking into %d placement subsets." %len(ptm))
            placement_tree_map.update(ptm)
            # self.root_problem.add_list_of_children(self.root_job_list)

            ''' For placement subsets create a placement subproblem, and decompose further'''
            for (p_key,p_tree) in ptm.items():
                assert isinstance(p_tree, PhylogeneticTree)
                # placement_problem  = SeppProblem(p_tree.leaf_node_names(), self.root_problem)
                placement_problem  = SeppProblem(p_tree.leaf_node_names(), part['problem'])
                placement_problem.subtree = p_tree
                placement_problem.label = "P_%s" %str(p_key)
                _LOG.debug("Placement subset %s has %d nodes" %(placement_problem.label,len(p_tree.leaf_node_names())))
                ''' Further decompose to alignment subsets '''
                alignment_tree_map = PhylogeneticTree(Tree(p_tree.den_tree)).decompose_tree(
                                            self.options.alignment_size,
                                            strategy=self.strategy,
                                            minSize = self.minsubsetsize,
                                            tree_map = {}, decomp_strategy = self.options.decomp_strategy)
                assert len(alignment_tree_map) > 0, ("Tree could not be decomposed"
                " given the following settings; strategy:%s minsubsetsize:%s alignmet_size:%s"
                %(self.strategy, self.minsubsetsize, self.options.alignment_size))

                _LOG.debug("Placement subset %s has %d alignment subsets: %s" %(placement_problem.label,len(alignment_tree_map.keys()),str(sorted(alignment_tree_map.keys()))))
                _LOG.debug("Placement subset %s has %d taxa:" %(placement_problem.label,sum([len(a_tree.leaf_node_names()) for a_tree in alignment_tree_map.values()])))
                for (a_key, a_tree) in alignment_tree_map.items():
                    assert isinstance(a_tree, PhylogeneticTree)
                    self.modify_tree(a_tree)
                    alignment_problem  = SeppProblem(a_tree.leaf_node_names(),
                                                      placement_problem)
                    alignment_problem.subtree = a_tree
                    alignment_problem.label = "A_%s_%s" %(str(p_key),str(a_key))

    def save_root_problem_to_pickle(self, filename=None):
        outf = open(os.path.join(storage_path,'root_problem.pkl'),'wb')
        pickle.dump(self.root_problem,outf)
        outf.close()

    def open_hmm_job_log(self):
        self.job_log = open(os.path.join(storage_path,'hmmbuild_job_log.txt'),'w')

    def close_hmm_job_log(self):
        self.job_log.close()

    def do_HMMbuild_jobs(self):
        '''
        Build all the HMMs that we will eventually need for this job and keep them organized.
        :return:
        '''
        self.master_hmm_reference={}
        job_log_line = '%(partition)s\t%(pp_name)s\t%(ap_name)s\t%(subaln_path)s\t%(hmm_path)s\n'

        for part in self.root_job_list.values():
            name = part['name']
            rp = part['problem']
            for pp in rp.children:

                for ap in pp.children:
                    ap_name = ap.label
                    fasta_fn = os.path.join(storage_path,ap_name + '_subaln.fas')
                    ap.write_subalignment_without_allgap_columns(fasta_fn)
                    hmm_fn = os.path.join(storage_path, ap_name + '.hmm')
                    bj = HMMBuildJob()
                    bj.setup(fasta_fn,hmm_fn)
                    args = {
                        'partition': name,
                        'pp_name': pp.label,
                        'ap_name': ap_name,
                        'subaln_path': fasta_fn,
                        'hmm_path': hmm_fn,
                        'job':bj
                    }
                    self.job_log.write(job_log_line % args)
                    self.master_hmm_reference[ap_name] = args.copy()

        for v in self.master_hmm_reference.values():
            JobPool().enqueue_job(v['job'])

        JobPool().wait_for_all_jobs()

    def build_hmms_and_save(self):
        _LOG.info('Preparing to initialize problems and save to folder: %s' % storage_path)
        self.build_subproblems_from_scratch()
        self.open_hmm_job_log()

        _LOG.info('Buildling HMM files')
        self.do_HMMbuild_jobs()

        _LOG.info('Saving root problem to pickle file to %s' % os.path.join(storage_path,'root_problem.pkl'))
        self.save_root_problem_to_pickle()

        self.close_hmm_job_log()

    def load_root_problem_from_pickled(self):
        ifile = open(os.path.join(storage_path,'root_problem.pkl'),'rb')
        self.root_problem = pickle.load(ifile)
        ifile.close()


