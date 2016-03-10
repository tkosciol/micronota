# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from os import makedirs
from os.path import join

from burrito.parameters import FlagParameter, ValuedParameter
from burrito.util import CommandLineApplication, ResultPath


class SignalP(CommandLineApplication):
    '''SignalP (version 4.1) application controller.'''
    _command = 'signalp'
    _valued_path_options = [
        # Logfile if -v is defined. Default: 'STDERR'
        '-l',
        # Specify temporary file directory. Default: /var/tmp
        '-T',
        # Make fasta file with mature sequence. Default: None
        '-m',
        # Make gff file of processed sequences. Default: 'Off'
        '-n'
    ]
    _valued_nonpath_options = [
        # Setting the output format ('short', 'long', 'summary' or 'all').
        # Default: 'short'
        '-f',
        # Graphics 'png' or 'png+eps'. Requires GNUPLOT. Default: 'Off'
        '-g',
        # Signal peptide networks to use ('best' or 'notm'). Default: 'best'
        '-s',
        # Organism type> (euk, gram+, gram-). Default: 'euk'
        '-t',
        # user defined D-cutoff for noTM networks
        '-u',
        # user defined D-cutoff for TM networks
        '-U',
        # Minimal predicted signal peptide length. Default: [10]
        '-M',
        # truncate to sequence length - 0 means no truncation. Default '70'
        '-c'
    ]
    _flag_options = [
        # Output this handy help message
        '-h',
        # Print SignalP version and exit
        '-V',
        # Verbose mode
        '-v',
        # Keep temporary directory. Default: 'Off'
        '-k',
        # web predictions. Default: 'Off'
        '-w'
    ]

    _parameters = {}
    _parameters.update({
        i: ValuedParameter(
            Prefix=i[0], Name=i[1:], Delimiter=' ',
            IsPath=True)
        for i in _valued_path_options})
    _parameters.update({
        i: ValuedParameter(
            Prefix=i[0], Name=i[1:], Delimiter=' ')
        for i in _valued_nonpath_options})
    _parameters.update({
        i: FlagParameter(
            Prefix=i[0], Name=i[1:])
        for i in _flag_options})
    _suppress_stderr = False

    def _accept_exit_status(self, exit_status):
        return exit_status == 0

    def _get_result_paths(self, data):
        result = {}
        # out_fp = re.sub(r"\>\ ", '', data[1])
        # result['output'] = ResultPath(out_fp, IsWritten=True)
        # result['output'] = sys.stdout

        # if `-k` flag is defined get temporaty file dir from `-T`
        if self.Parameters['-k'].isOn():
            if self.Parameters['-T'].isOn():
                tmp_fp = self._absolute(self.Parameters['-T'].Value)
            else:
                # taken from default definition in `_valued_path_options`
                tmp_fp = '/var/tmp/'
            result['tmp'] = ResultPath(Path=tmp_fp, IsWritten=True)

        # get log, fasta and gff file paths
        for option, vals in [('-l', 'log'), ('-m', 'fasta'), ('-n', 'gff')]:
            if self.Parameters[option].isOn():
                fp = self._absolute(self.Parameters[option].Value)
                result[vals] = ResultPath(Path=fp, IsWritten=True)

        # # SKIP THIS FOR NOW
        # # get png (and eps) files
        # if self.Parameters['-g'].isOn():
        #     # get inp_fp GI_IDs
        #     # # HERE
        #     gis = []
        #     # get png files
        #     for gi in gis:
        #         result['png']
        #     if self.Parameters['-g'].Value is 'gff+eps':
        #         # get eps files
        #         result['eps']

        return result


def predict_signal(in_fp, out_dir, prefix, params=None):
    '''Predict signal peptide cleavage sites for the input file.

    Notes
    -----
    It will create an output file, depending on the selected parameter:
        A. short
        B. long
        C. summary
        D. all
    SignalP accepts any input and does not raise input file errors.
    Please check for correct input before running SignalP.

    Parameters
    ----------
    in_fp : str
        input file path
    out_dir : str
        output file directory path
    prefix : str
        name of the output file
    params : dict
        Other command line parameters for SignalP. key is the option
        (e.g. "-t") and value is the value for the option (e.g. "euk").
        If the option is a flag, set the value to None.

    Returns
    -------
    burrito.util.CommandLineAppResult
        It contains opened file handlers of stdout, stderr, and the
        output files, which can be accessed in a dict style with the
        keys of "StdOut", "StdErr", "output", "tmp" (if specified),
        "gff" (if specified), "fasta" (if specified), "png" & "eps" (if spec.)
        and "log" (if specified). The exit status of the run can be similarly
        fetched with the key "ExitStatus".
    '''
    # create dir if does not exist
    makedirs(out_dir, exist_ok=True)

    if params is None:
        params = {}

    # determine suffix on the basis of '-f' option
    if '-f' in params:
        out_suffix = "_".join(["sp", params['-f']])
    else:
        out_suffix = "sp_short"

    # change -m, -n and -l options to include path if there isn't one!
    for i in ['-l', '-m', '-n']:
        if i in params:
            if '/' not in params[i]:
                params[i] = join(out_dir, params[i])

    out_fp = join(out_dir, '.'.join([prefix, out_suffix]))

    app = SignalP(InputHandler='_input_as_paths', params=params)
    return app([in_fp, out_fp])
