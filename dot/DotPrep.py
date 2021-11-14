import os
import warnings
from tempfile import TemporaryDirectory

from Bio import BiopythonParserWarning
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature

from .dot_prep_orig import run as dotprep_run
from .Nucmer import Nucmer


class DotPrepArgs:
    def __init__(self, delta, output_filename, unique_length=10000, max_overview_alignments=1000):
        self.delta = delta
        self.unique_length = unique_length
        self.out = output_filename
        self.overview = max_overview_alignments


class DotPrep:
    def __init__(self):
        self.nucmer = Nucmer()

    def run_shell(self, outdir: str, fasta_ref: str, fasta_qry: str, mincluster: int = 65, gbk_ref: str = None, gbk_qry: str = None) -> str:
        """
        Run nucmer on FASTAs, run DotPrep, and optionally create Dot annotation files based from GenBank files.

        :param outdir: path to output directory
        :param fasta_ref: path to first assembly FASTA
        :param fasta_qry: path to second assembly FASTA
        :param gbk_ref: path to first GenBank file (optional)
        :param gbk_qry: path to second GenBank file (optional)
        :param mincluster: sets the minimum length of a cluster of matches
        """
        assert not os.path.isdir(outdir), f'{outdir=} already exists!'
        os.makedirs(outdir)

        coords, index = self._run(workdir=outdir, fasta_ref=fasta_ref, fasta_qry=fasta_qry, mincluster=mincluster)

        # load reference GenBank
        if gbk_ref:
            with open(f'{outdir}/out.ref.annotations', 'w') as f:
                f.write(self.gbk_to_annotation_file(gbk=gbk_ref, is_ref=True))

        # load query GenBank
        if gbk_qry:
            with open(f'{outdir}/out.qry.annotations', 'w') as f:
                f.write(self.gbk_to_annotation_file(gbk=gbk_qry, is_ref=False, index=index))

        return f'Complete success. {outdir=}'

    def run_python(self, fasta_ref: str, fasta_qry: str, mincluster: int = 65) -> (str, str):
        """
        Returns coords and index files as string.

        :param fasta_ref: path to first assembly FASTA
        :param fasta_qry: path to second assembly FASTA
        :param mincluster: sets the minimum length of a cluster of matches
        """
        with TemporaryDirectory() as workdir:
            coords, index = self._run(workdir=workdir, fasta_ref=fasta_ref, fasta_qry=fasta_qry, mincluster=mincluster)

        return coords, index

    def _run(self, workdir: str, fasta_ref: str, fasta_qry: str, mincluster: int) -> (str, str):
        for fasta in [fasta_ref, fasta_qry]:
            assert os.path.isfile(fasta), f'path is invalid: {fasta}'

        # run nucmer
        delta_path = self.nucmer.align(fasta_ref=fasta_ref, fasta_qry=fasta_qry, workdir=workdir, arguments={'--mincluster': mincluster})

        # run dotprep
        dotprep_run(DotPrepArgs(delta=delta_path, output_filename=f'{workdir}/out'))

        # verify results
        with open(f'{workdir}/out.coords') as f:
            coords = f.read()
        with open(f'{workdir}/out.coords.idx') as f:
            index = f.read()
        assert coords.split('\n', 1)[0] == 'ref_start,ref_end,query_start,query_end,ref', f'DotPrep failed: coords incorrect. {coords=}'
        assert index.split('\n', 1)[0] == '#ref', f'DotPrep failed: index incorrect. {index=}'

        return coords, index

    @classmethod
    def gbk_to_annotation_file(cls, gbk: str, is_ref: bool, index: str = None) -> str:
        """
        Turn gbk into Dot-compliant file format.

        :param gbk: path to GenBank file that corresponds to fasta_ref
        :param is_ref: are the annotations for the reference- or the query sequence?
        :param index: only required if is_ref=False; path to out.index or string
        :return: dot-compliant annotations string
        """
        if is_ref:
            result = 'ref,ref_start,ref_end,name,strand\n'
        else:
            result = 'query,query_start,query_end,name,strand\n'

        for annotation in cls.gbk_to_annotation_dict(gbk=gbk, is_ref=is_ref, index=index):
            scf_id, start, end, locus_tag, strand = annotation.values()
            result += f"{scf_id},{start},{end},{locus_tag},{strand}\n"

        return result

    @classmethod
    def gbk_to_annotation_dict(cls, gbk: str, is_ref: bool, index: str = None) -> [dict]:
        """
        Turn gbk into Dot-compliant dict format.

        :param gbk: path to GenBank file that corresponds to fasta_ref
        :param is_ref: are the annotations for the reference- or the query sequence?
        :param index: only required if is_ref=False; path to out.index or string
        :return: dot-compliant annotations string
        """
        if is_ref:
            return cls.gbk_to_annotation_ref(gbk=gbk)
        else:
            assert index is not None, f'With is_ref=False, index must be specified. {index=}'
            if '\n' not in index:
                with open(index) as f:
                    index = f.read()

            flipped_scaffolds = cls._get_flipped_scaffolds(index)
            return cls.gbk_to_annotation_qry(gbk=gbk, flipped_scaffolds=flipped_scaffolds)

    @staticmethod
    def gbk_to_annotation_qry(gbk: str, flipped_scaffolds: {str}) -> [dict]:
        """
        Turn gbk into Dot-compliant dict format (query).

        :param gbk: path to genbank file that corresponds to fasta_ref
        :param flipped_scaffolds: dictionary of flipped scaffolds. only supply if is_ref=False
        :return: dot-compliant annotations list of dictionaries
        """
        assert type(flipped_scaffolds) is dict, f'To load query annotations, supply set of flipped scaffolds! {type(flipped_scaffolds)}'

        result = []
        with warnings.catch_warnings(), open(gbk) as input_handle:
            warnings.simplefilter('ignore', BiopythonParserWarning)

            for scf in SeqIO.parse(input_handle, "genbank"):

                if scf.id not in flipped_scaffolds:
                    # skip scaffolds that are not in the dot plot
                    continue

                len_scf = len(scf)
                is_flipped = flipped_scaffolds[scf.id] == '-'

                loci = set()  # avoid duplicate genes

                for f in scf.features:
                    if 'locus_tag' not in f.qualifiers:
                        continue

                    start, end = DotPrep.__get_position(f=f)

                    if (start, end) in loci:
                        continue

                    if is_flipped:
                        # mirror the genes on this scaffold
                        start, end = len_scf - end, len_scf - start

                    result.append({
                        'query': scf.id,
                        'query_start': start,
                        'query_end': end,
                        'name': f.qualifiers['locus_tag'][0],
                        'strand': '+' if f.strand > 0 else '-'
                    })

                    loci.add((start, end))

        return result

    @staticmethod
    def gbk_to_annotation_ref(gbk: str) -> [dict]:
        """
        Turn gbk into Dot-compliant dict format (reference).

        :param gbk: path to genbank file that corresponds to fasta_ref
        :return: dot-compliant annotations list of dictionaries
        """

        result = []
        with warnings.catch_warnings(), open(gbk) as input_handle:
            warnings.simplefilter('ignore', BiopythonParserWarning)

            for scf in SeqIO.parse(input_handle, "genbank"):
                loci = set()  # {(1, 300), (444, 600), ...}
                for f in scf.features:
                    if 'locus_tag' not in f.qualifiers:
                        continue

                    start, end = DotPrep.__get_position(f=f)

                    if (start, end) in loci:
                        continue

                    result.append({
                        'ref': scf.id,
                        'ref_start': start,
                        'ref_end': end,
                        'name': f.qualifiers['locus_tag'][0],
                        'strand': '+' if f.strand > 0 else '-'
                    })

                    loci.add((start, end))

        return result

    @staticmethod
    def _get_flipped_scaffolds(index: str) -> {str: str}:
        query_scaffolds = index.split('matching_refs\n', 1)[1].split('#overview', 1)[0].rstrip()

        if query_scaffolds == '':
            return {}  # no alignments

        def extract_(line):
            scf_name, _, orientation, _ = line.split(',', 3)
            return scf_name, orientation

        return dict(extract_(line) for line in query_scaffolds.split('\n'))

    @staticmethod
    def __get_position(f: SeqFeature) -> (int, int):
        location = f.location if f.location_operator != 'join' else f.location.parts[0]
        return (location.nofuzzy_start, location.nofuzzy_end)


def main():
    import fire

    dotprep = DotPrep()

    fire.Fire({
        'run': dotprep.run_shell,
        'create_annotations': dotprep.gbk_to_annotation_file,
    })


if __name__ == '__main__':
    main()
