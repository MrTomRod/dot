import os
import shutil
from unittest import TestCase
from shutil import rmtree
from dot import dot_prep_orig, DotPrep, DotPrepArgs, Nucmer

# fasta_ref = 'tests/fna/FAM19036.fna'
# gbk_ref = 'tests/gbk/FAM19036.gbk'
# fasta_qry = 'tests/fna/FAM14217.fna'
# gbk_qry = 'tests/gbk/FAM14217.gbk'

fasta_ref = 'tests/fna/FAM19036.fna'
gbk_ref = 'tests/gbk/FAM19036.gbk'
fasta_qry = 'tests/fna/FAM19036i.fna'
gbk_qry = 'tests/gbk/FAM19036i.gbk'

# fasta_ref = 'tests/fakedelta/REFERENCE.fna'
# gbk_ref = 'tests/fakedelta/REFERENCE.gbk'
# fasta_qry = 'tests/fakedelta/QUERY.fna'
# gbk_qry = 'tests/fakedelta/QUERY.gbk'

nucmer = Nucmer()
dotprep = DotPrep()


def rm_dir(dir: str) -> str:
    if os.path.isdir(dir):
        shutil.rmtree(dir)
    return dir


class Test(TestCase):
    def test_nucmer(self):
        tmp_dir = '/tmp/test_nucmer'
        if os.path.isdir(tmp_dir):
            rmtree(tmp_dir)
        os.makedirs(tmp_dir)
        res = nucmer.align(
            fasta_ref=fasta_ref, fasta_qry=fasta_qry,
            workdir=tmp_dir,
            arguments={'--mincluster': 100}
        )
        assert os.path.isfile(res)

    def test_get_flipped_scaffolds(self):
        index = '''#ref
ref,ref_length,matching_queries
FAM19036-p1-1_scf0001,3622795,FAM19036-i1-1_scf21~FAM19036-i1-1_scf27~FAM19036-i1-1_scf10~FAM19036-i1-1_scf16~FAM19036-i1-1_scf24~FAM19036-i1-1_scf15~FAM19036-i1-1_scf5~FAM19036-i1-1_scf23~FAM19036-i1-1_scf8~FAM19036-i1-1_scf1~FAM19036-i1-1_scf4~FAM19036-i1-1_scf25~FAM19036-i1-1_scf29~FAM19036-i1-1_scf11~FAM19036-i1-1_scf12~FAM19036-i1-1_scf2~FAM19036-i1-1_scf30~FAM19036-i1-1_scf37~FAM19036-i1-1_scf32~FAM19036-i1-1_scf6~FAM19036-i1-1_scf22~FAM19036-i1-1_scf17~FAM19036-i1-1_scf7~FAM19036-i1-1_scf26~FAM19036-i1-1_scf20~FAM19036-i1-1_scf35~FAM19036-i1-1_scf33~FAM19036-i1-1_scf18~FAM19036-i1-1_scf3~FAM19036-i1-1_scf19~FAM19036-i1-1_scf28~FAM19036-i1-1_scf13~FAM19036-i1-1_scf14~FAM19036-i1-1_scf31~FAM19036-i1-1_scf9
#query
query,query_length,orientation,bytePosition_unique,bytePosition_repetitive,bytePosition_end,unique_matching_refs,matching_refs
FAM19036-i1-1_scf34,12253,+,4477,28,290,,FAM19036-p1-1_scf0001
FAM19036-i1-1_scf36,6509,+,9609,28,242,,FAM19036-p1-1_scf0001
FAM19036-i1-1_scf38,2147,+,9911,28,196,,FAM19036-p1-1_scf0001
FAM19036-i1-1_scf4,277130,-,1071,165,361,FAM19036-p1-1_scf0001,FAM19036-p1-1_scf0001
FAM19036-i1-1_scf2,352264,-,4877,172,809,FAM19036-p1-1_scf0001,FAM19036-p1-1_scf0001
FAM19036-i1-1_scf18,63098,+,2873,72,238,FAM19036-p1-1_scf0001,FAM19036-p1-1_scf0001
#overview
ref_start,ref_end,query_start,query_end,ref,query,tag
269129,475154,0,205999,FAM19036-p1-1_scf0001,FAM19036-i1-1_scf2,unique
        '''
        flipped_scaffolds = dotprep._get_flipped_scaffolds(index=index)
        self.assertEqual(
            flipped_scaffolds,
            {'FAM19036-i1-1_scf34': '+', 'FAM19036-i1-1_scf36': '+', 'FAM19036-i1-1_scf38': '+',
             'FAM19036-i1-1_scf4': '-', 'FAM19036-i1-1_scf2': '-', 'FAM19036-i1-1_scf18': '+'}
        )

    def test_get_flipped_scaffolds_empty(self):
        index = '''#ref
#query
query,query_length,orientation,bytePosition_unique,bytePosition_repetitive,bytePosition_end,unique_matching_refs,matching_refs
#overview
ref_start,ref_end,query_start,query_end,ref,query,tag
        '''
        flipped_scaffolds = dotprep._get_flipped_scaffolds(index=index)
        self.assertEqual(flipped_scaffolds, {})

    def test_dotprep_py(self):
        coords, index = dotprep.run_python(fasta_ref=fasta_ref, fasta_qry=fasta_qry, mincluster=65)

    def test_dotprep_shell(self):
        msg = dotprep.run_shell(
            outdir=rm_dir('tests/outdir'),
            fasta_ref=fasta_ref, fasta_qry=fasta_qry,
            gbk_ref=gbk_ref, gbk_qry=gbk_qry,
            mincluster=65
        )
        print(msg)

    def test_getqrc(self):
        header_lines_by_query, lines_by_query = dot_prep_orig.getQueryRefCombinations('tests/fakedelta/out.delta')
        print(header_lines_by_query, lines_by_query)

    def test_gbk_to_annotation_file(self):
        annotations = DotPrep.gbk_to_annotation_file(gbk_ref, is_ref=True)
        print(annotations[:200])

    def test_gbk_to_annotation_qry(self):
        annotations = DotPrep.gbk_to_annotation_file(
            gbk='tests/gbk/FAM19036i.gbk',
            is_ref=False,
            index='tests/q_FAM19036i_r_FAM19036/output.coords.idx',
        )
        print(annotations.split('\n', 4)[:-1])
        with open('tests/q_FAM19036i_r_FAM19036/output.annotations', 'w') as f:
            f.write(annotations)

    def test_gbk_to_annotation_ref(self):
        annotations = DotPrep.gbk_to_annotation_file(
            gbk='tests/gbk/FAM19036.gbk',
            is_ref=True
        )
        print(annotations.split('\n', 4)[:-1])
        with open('tests/q_FAM19036i_r_FAM19036/output.annotations', 'w') as f:
            f.write(annotations)

    def test_dotprep_delta(self):
        dot_prep_orig.run(DotPrepArgs(delta='tests/fakedelta/out.delta.original', output_filename='tests/fakedelta/out'))

    def test_dotprep_synth_a(self):
        Nucmer().align(fasta_ref='tests/synth_data/REFERENCE.fna', fasta_qry='tests/synth_data/QUERY.fna', workdir='tests/synth_data')

    def test_dotprep_synth(self):
        dot_prep_orig.run(DotPrepArgs(delta='tests/synth_data/out.delta', output_filename='tests/synth_data/out'))
