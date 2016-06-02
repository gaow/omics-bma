#! /usr/bin/env python3
__author__ = "Gao Wang"
__copyright__ = "Copyright 2016, Stephens lab"
__email__ = "gaow@uchicago.edu"
__license__ = "MIT"
__version__ = "0.1.0"
import sys, os, subprocess, shlex, re, \
    datetime, gzip, time, bz2, yaml, collections
from io import StringIO
from contextlib import contextmanager
import pandas as pd
from pysos.utils import RuntimeEnvironments

class Environment(RuntimeEnvironments):
    _instance = None
    def __new__(cls, *args, **kwargs):
        if not cls._instance:
            cls._instance = super(Environment, cls).__new__(cls)
        return cls._instance

    def __init__(self):
        super().__init__()
        self.term_width = None
        self.__width_cache = 1
        self.path = {'PATH':"{}:{}".format(os.getcwd(), os.environ["PATH"])}
        self.debug = False
        self.quiet = False
        self.colors = "#377EB8 #E41A1C #4DAF4A #984EA3 #FFD92F #FF7F00 #F781BF " \
                      "#8DD3C7 #B3B3B3 #000000 #56B4E9 #BC80BD #FDB462 #350E20 #8A9045 #800000".split()

    def error(self, msg = None, show_help = False, exit = False):
        if msg is None:
            sys.stderr.write('\n')
            return
        if type(msg) is list:
            msg = ' '.join(map(str, msg))
        else:
            msg = str(msg)
        start = '\n' if msg.startswith('\n') else ''
        end = '\n' if msg.endswith('\n') else ''
        msg = msg.strip()
        if exit:
            tag = 'ERROR'
        else:
            tag = 'WARNING'
        sys.stderr.write(start + "\033[1;40;33m{}: {}\033[0m\n".format(tag, msg) + end)
        if show_help:
            self.log("Type '--help' for help message")
            sys.exit()
        if exit:
            sys.exit()

    def log(self, msg = None, flush=False):
        if self.debug or self.quiet:
            return
        if msg is None:
            sys.stderr.write('\n')
            return
        if type(msg) is list:
            msg = ' '.join(map(str, msg))
        else:
            msg = str(msg)
        start = "{0:{width}}".format('\r', width = self.__width_cache + 10) + "\r" if flush else ''
        end = '' if flush else '\n'
        start = '\n' + start if msg.startswith('\n') else start
        end = end + '\n' if msg.endswith('\n') else end
        msg = msg.strip()
        sys.stderr.write(start + "\033[1;40;32mMESSAGE: {}\033[0m".format(msg) + end)
        self.__width_cache = len(msg)

env = Environment()

def openFile(filename):
    if filename.lower().endswith('.tar.gz') or filename.lower().endswith('.tgz'):
        raise RuntimeError('Please decompress {} before reading.'.format(filename))
    if filename.lower().endswith('.gz'):
        return gzip.open(filename, 'rb')
    elif filename.lower().endswith('.bz2'):
        return bz2.BZ2File(filename, 'rb')
    else:
        # text file
        # because readline() from gzip.open will be byte, not string, we should return
        # binary here in order to process them equally in order for things to work
        # correctly under python 3
        return open(filename, 'rb')

def rename_tmp(filename):
    '''Temporary output of filename'''
    # turn path/filename.ext to path/filename_tmp???.ext, where ??? is
    # the process ID to avoid two processes writing to the same temp
    # files. That is to say, if two processes are working on the same step
    # they will produce different temp files, and the final results should
    # still be valid.
    return '_tmp{}.'.format(os.getpid()).join(filename.rsplit('.', 1))

def rename_stamp(filename):
    '''Stamped output of filename'''
    # turn path/filename.ext to path/filename-???.ext where ???
    # is the date and if available the git branch name
    # Figure out the date
    time = datetime.datetime.now()
    date = time.strftime("%Y") + time.strftime("%m") + time.strftime("%d")
    try:
        branch = runCommand('git rev-parse --abbrev-ref HEAD')[0]
    except:
        branch = None
    stamp = '_{}{}.'.format(date, '_{}'.format(branch) if branch is not None else '')
    return stamp.join(filename.rsplit('.', 1))

def get_value_type(x):
    if x is None:
        return 'null'
    if x.lower() in ['na','nan','null','none']:
        x = 'nan'
    try:
        int(x)
        return 'int'
    except:
        try:
            float(x)
            return 'float'
        except:
            return 'string'

def is_null(x):
    if x is None:
        return True
    if type(x) is str:
        if x.lower() in ['na','nan','null','none']:
            return True
    return False

def str2num(x):
    if type(x) is str:
        try:
            return int(x)
        except ValueError:
            try:
                return float(x)
            except ValueError:
                return x
    else:
        return x

class StdoutCapturer(list):
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self._stringio = StringIO()
        return self

    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        sys.stdout = self._stdout

@contextmanager
def stdoutRedirect(to=os.devnull):
    '''
    import os

    with stdoutRedirect(to=filename):
        print("from Python")
        os.system("echo non-Python applications are also supported")
    '''
    fd = sys.stdout.fileno()

    ##### assert that Python and C stdio write using the same file descriptor
    ####assert libc.fileno(ctypes.c_void_p.in_dll(libc, "stdout")) == fd == 1

    def _redirect_stdout(to):
        sys.stdout.close() # + implicit flush()
        os.dup2(to.fileno(), fd) # fd writes to 'to' file
        sys.stdout = os.fdopen(fd, 'w') # Python writes to fd

    with os.fdopen(os.dup(fd), 'w') as old_stdout:
        with open(to, 'a') as file:
            _redirect_stdout(to=file)
        try:
            yield # allow code to be run with the redirected stdout
        finally:
            _redirect_stdout(to=old_stdout) # restore stdout.
                                            # buffering and flags such as
                                            # CLOEXEC may be different

def runCommand(cmd, instream = None, msg = '', upon_succ = None, show_stderr = False, return_zero = True):
    if isinstance(cmd, str):
        cmd = shlex.split(cmd)
    popen_env = os.environ.copy()
    popen_env.update(env.path)
    try:
        tc = subprocess.Popen(cmd, stdin = subprocess.PIPE,
                              stdout = subprocess.PIPE, stderr = subprocess.PIPE,
                              env=popen_env)
        if instream:
            instream = instream.encode(sys.getdefaultencoding())
            out, error = tc.communicate(instream)
        else:
            out, error = tc.communicate()
        out = out.decode(sys.getdefaultencoding())
        error = error.decode(sys.getdefaultencoding())
        if return_zero:
            if tc.returncode < 0:
                raise ValueError ("Command '{0}' was terminated by signal {1}".format(cmd, -tc.returncode))
            elif tc.returncode > 0:
                raise ValueError ("{0}".format(error))
        if error.strip() and show_stderr:
            env.logger.error(error)
    except OSError as e:
        raise OSError ("Execution of command '{0}' failed: {1}".format(cmd, e))
    # everything is OK
    if upon_succ:
        # call the function (upon_succ) using others as parameters.
        upon_succ[0](*(upon_succ[1:]))
    return out.strip(), error.strip()

def is_empty(v):
    if type(v) == pd.DataFrame:
        return v.empty
    else:
        if v:
            return False
        else:
            return True

def replace_non_alnum(x):
    trans = ''.join(chr(c) if chr(c).isalnum() else '_' for c in range(256))
    x.translate(trans)

class Timer(object):
    def __init__(self, verbose=False):
        self.verbose = verbose

    def __enter__(self):
        self.start = time.time()
        return self

    def __exit__(self, *args):
        self.end = time.time()
        self.secs = self.end - self.start
        self.msecs = self.secs * 1000  # millisecs
        if self.verbose:
            print('elapsed time: %.03f ms' % self.msecs)

def flatten_dict(d):
    items = []
    for k, v in d.items():
        if isinstance(v, collections.MutableMapping):
            items.extend(flatten_dict(v).items())
        else:
            items.append((k, v))
    return dict(items)

def str2list(value):
    if value is None:
        return []
    else:
        return [x.strip() for x in re.split(" |\+|,", value) if x.strip()]

def strip_dict(data, mapping = dict, into_list = False):
    if not isinstance(data, mapping):
        return data
    mapping_null = mapping()
    new_data = mapping()
    for k, v in data.items():
        if isinstance(v, collections.Mapping):
            v = strip_dict(v, mapping, into_list)
        if isinstance(v, list) and into_list:
            v = [strip_dict(x, mapping, into_list) for x in v]
        if not is_null(v) and v != mapping_null:
            new_data[k] = v
    return new_data

def dict2str(value, replace = []):
    out = StringIO()
    yaml.dump(strip_dict(value, into_list = True), out, default_flow_style=False)
    res = out.getvalue()
    out.close()
    for item in replace:
        res = res.replace(item[0], item[1])
    return res

def test_almost_equal_recursive(x, y, level = 6, level_cutoff = 3):
    try:
      np.testing.assert_almost_equal(x, y, level)
    except AssertionError:
      if level <= level_cutoff:
          raise
      else:
          level = level - 1
          test_almost_equal_recursive(x, y, level)
    return level
