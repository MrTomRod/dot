import os
from subprocess import call, run, PIPE
from .utils import is_installed


class Nucmer:
    """
    Wrapper for nucmer3: http://mummer.sourceforge.net/
    Requires nucmer to be installed.
    """

    def __init__(self, nucmer_path='nucmer'):
        self.nucmer_path = nucmer_path
        assert is_installed(self.nucmer_path), f'Nucmer is not installed! {self.nucmer_path}'

        # nucmer works
        assert call([self.nucmer_path, '--version'], stderr=PIPE, stdout=PIPE) == 0, 'Nucmer is not executable!'

    def align(self, fasta_ref: str, fasta_qry: str, workdir: str, arguments: dict = None, clean_input: bool = True) -> str:
        """
        :param fasta_ref: path to first assembly fasta
        :param fasta_qry: path to second assembly fasta
        :param arguments: dict with additional arguments for nucmer, e.g. {'--mincluster': 100}
        :param clean_input: if true, remove prefixes such as 'gnl|Prokka|' or 'gnl|extdb|' from input fnas
        :return: str: path to out.delta
        """
        assert os.path.isdir(workdir), f'work_dir does not exist: {workdir}'
        result_path = f'{workdir}/out.delta'
        assert not os.path.isfile(result_path), f'result_path already exists: {result_path}'

        if clean_input:
            # create clean input fastas in temp dir
            fasta_ref = self.remove_fasta_prefix(in_fasta=fasta_ref, out_fasta=f'{workdir}/ref.fasta.cleaned')
            fasta_qry = self.remove_fasta_prefix(in_fasta=fasta_qry, out_fasta=f'{workdir}/qry.fasta.cleaned')

        cmd = [self.nucmer_path]
        if arguments is not None:
            for arg, val in arguments.items():
                cmd.extend([str(arg), str(val)])
        cmd.extend([os.path.abspath(fasta_ref), os.path.abspath(fasta_qry)])

        subprocess = run(
            cmd,
            cwd=workdir,
            stdout=PIPE, stderr=PIPE, encoding='ascii')

        error_message = f'Nucmer occurred with fasta_ref={fasta_ref} and fasta_qry={fasta_qry}; stdout={subprocess.stdout}; stderr={subprocess.stderr}'

        assert subprocess.returncode == 0, error_message
        assert os.path.isfile(result_path), error_message
        return result_path

    @staticmethod
    def remove_fasta_prefix(in_fasta: str, out_fasta: str) -> str:
        """
        Remove prefixes such as 'gnl|Prokka|' or 'gnl|extdb|' from protein/contig identifiers.
        Also removes comments that follow after the first blank.

        E.g. '>ncbi|prefix|gene_0001 comment' ➜ '>gene_0001'

        :param in_fasta: path to input fasta
        :param out_fasta: path to output fasta
        :type out_fasta: same as out_fasta
        """
        assert os.path.isfile(in_fasta), f'file does not exist: {in_fasta}'

        with open(in_fasta) as in_f, open(out_fasta, 'w') as out_f:
            for line in in_f:
                out_f.write(
                    '>' + line[1:-1].split(' ', 1)[0].rsplit('|', 1)[-1] + '\n'
                    if line.startswith('>')
                    else line
                )
        return out_fasta
