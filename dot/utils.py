import os


def is_installed(program):
    """
    Test if a program is installed.

    :param program: path to executable or command
    :return: if program executable: program; else None
    """

    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:  # check if path to program is valid
        return is_exe(program)
    else:  # check if program is in PATH
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
        return False
