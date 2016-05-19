#! /usr/bin/env python3
__author__ = "Gao Wang"
__copyright__ = "Copyright 2016, Stephens lab"
__email__ = "gaow@uchicago.edu"
__license__ = "MIT"
__version__ = "0.1.0"
import sys, os, subprocess, shutil, glob, shlex, re, \
    hashlib, datetime, gzip, time, bz2, array, \
    tarfile, platform
from distutils.dir_util import mkpath
from io import StringIO
from contextlib import contextmanager
import itertools
from collections import OrderedDict, defaultdict, Counter, MutableMapping
import pandas as pd

class Environment:
    def __init__(self):
        self.term_width = None
        self.__width_cache = 1
        self.path = {'PATH':"{}:{}".format(os.getcwd(), os.environ["PATH"])}
        self.debug = False
        self.quiet = False
        self.colors = "#377EB8 #E41A1C #4DAF4A #984EA3 #FFD92F #FF7F00 #F781BF #8DD3C7 #B3B3B3 #000000 #56B4E9 #BC80BD #FDB462 #350E20 #8A9045 #800000".split()

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

class PrettyPrinter:
    def __init__(self, delimiter=None, max_width={}, cache_size=200):
        ''' delimiter: use specified field to separate fields
            max_width: a dictionary of {col: max_width} to change long
                text to START...END
            cache: only use the first cache lines to get column width
        '''
        self.width = []
        self.rows = []
        self.max_width = max_width
        self.cache_size = cache_size
        # if a delimiter is specified, use it
        if delimiter is not None:
            self.delimiter = delimiter.replace(r'\t', '\t')
            self.write = self.direct_print
            self.write_rest = self.direct_print_rest
        elif max_width:
            self.delimiter = '\t'
            self.write = self.cached_trim_print
            self.write_rest = self.cached_trim_print_rest
        else:
            self.delimiter = '\t'
            self.write = self.cached_print
            self.write_rest = self.cached_print_rest

    #
    # MODE 1: direct print
    #
    def direct_print(self, data):
        '''print data directly using specified delimiter'''
        print(self.delimiter.join([x for x in data]))

    def direct_print_rest(self):
        '''No cache so do nothing'''
        pass

    #
    # MODE 2: cached, trimmed print
    #
    def cached_trim_print(self, data):
        '''Use cache, figure out column width'''
        trimmed = {}
        for c,m in list(self.max_width.items()):
            if len(data[c]) > m:
                trimmed[c] = data[c][: m // 3] + '...' + data[c][ - (m - m // 3 - 3):]
        if trimmed:
            trimmed_data = [x for x in data]
            for c,txt in list(trimmed.items()):
                trimmed_data[c] = txt
        else:
            trimmed_data = data
        #
        self.rows.append(trimmed_data)
        if not self.width:
            self.width = [len(x) for x in trimmed_data]
        else:
            self.width = [max(y, len(x)) for y,x in zip(self.width, trimmed_data)]
        # cache size exceeds, use collected width and stop checking
        if len(self.rows) > self.cache_size:
            self.cached_trim_print_rest()
            # change print mode
            self.write = self.uncached_trim_print

    def cached_trim_print_rest(self):
        '''Print and clear cache'''
        if not self.rows:
            return
        # do not ljust the last column. This avoids unnecessary spaces
        # at the end of each line
        self.width[-1] = 0
        # print everything in cache
        print('\n'.join([
            self.delimiter.join(
                [col.ljust(width) for col, width in zip(row, self.width)])
            for row in self.rows]))
        # clear cache
        self.rows = []

    def uncached_trim_print(self, data):
        trimmed = {}
        for c,m in list(self.max_width.items()):
            if len(data[c]) > m:
                trimmed[c] = data[c][: m // 3] + '...' + data[c][ - (m - m // 3 - 3):]
        if trimmed:
            trimmed_data = [x for x in data]
            for c,txt in list(trimmed.items()):
                trimmed_data[c] = txt
        else:
            trimmed_data = data
        #
        print(self.delimiter.join(
            [col.ljust(width) for col, width in zip(trimmed_data, self.width)]))

    #
    # MODE 3: cached, untrimmed print
    #
    def cached_print(self, data):
        self.rows.append(data)
        if not self.width:
            self.width = [len(x) for x in data]
        else:
            self.width = [max(y, len(x)) for y,x in zip(self.width, data)]
        # cache size exceeds, use collected width and stop checking
        if len(self.rows) > self.cache_size:
            self.cached_print_rest()
            # change print mode
            self.write = self.uncached_trim_print

    def cached_print_rest(self):
        if not self.rows:
            return
        # do not ljust the last column. This avoids unnecessary spaces
        # at the end of each line
        self.width[-1] = 0
        print('\n'.join([
            self.delimiter.join(
                [col.ljust(width) for col, width in zip(row, self.width)])
            for row in self.rows]))
        self.rows = []

    def uncached_print(self, data):
        print(self.delimiter.join(
            [col.ljust(width) for col, width in zip(data, self.width)]))

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

def calculateMD5(filename, partial=False):
    filesize = os.path.getsize(filename)
    # calculate md5 for specified file
    md5 = hashlib.md5()
    block_size = 2**20  # buffer of 1M
    try:
        if (not partial) or filesize < 2**26:
            with open(filename, 'rb') as f:
                while True:
                    data = f.read(block_size)
                    if not data:
                        break
                    md5.update(data)
        else:
            count = 64
            # otherwise, use the first and last 500M
            with open(filename, 'rb') as f:
                while True:
                    data = f.read(block_size)
                    count -= 1
                    if count == 32:
                        f.seek(-2**25, 2)
                    if not data or count == 0:
                        break
                    md5.update(data)
    except IOError as e:
        sys.exit('Failed to read {}: {}'.format(filename, e))
    return md5.hexdigest()

def compressFile(infile, outfile):
    '''Compress a file from infile to outfile'''
    with open(infile, 'rb') as input, gzip.open(outfile, 'wb') as output:
            buffer = input.read(100000)
            while buffer:
                output.write(buffer)
                buffer = input.read(100000)
    return outfile

def decompressGzFile(filename, inplace=True, force=False, md5=None):
    '''Decompress a file.gz and return file if needed'''
    if filename.lower().endswith('.tar.gz') or filename.lower().endswith('.tar.bz2'):
        dest_files = []
        mode = 'r:gz'
        with tarfile.open(filename, mode) as tar:
            # only extract files
            path = os.path.dirname(filename)
            files = [x.name for x in tar.getmembers() if x.isfile()]
            for f in files:
                dest_file = os.path.join(path, os.path.basename(f))
                dest_files.append(dest_file)
                if not os.path.isfile(dest_file):
                    tar.extract(f, path)
        return dest_files
    elif filename.lower().endswith('.gz'):
        new_filename = filename[:-3]
        # if the decompressed file exists, and is newer than the .gz file, ignore
        if os.path.isfile(new_filename) and not force:
            if md5 is not None and md5 != calculateMD5(new_filename, partial=True):
                env.error('MD5 signature mismatch: {} (signature: {} calculated: {})'
                    .format(new_filename, md5, calculateMD5(new_filename, partial=True)))
            else:
                env.log('Reusing existing decompressed file {}'.format(new_filename))
                return new_filename
        #
        # check if dest_dir is writable
        dest_dir = os.path.dirname(filename)
        if not os.access(dest_dir, os.W_OK):
            # if we are decompressing files from a read-only shared repository
            # write to local_resource
            if os.path.abspath(dest_dir).startswith(os.path.abspath(env.shared_resource)):
                new_filename = '{}/{}'.format(env.local_resource,
                    os.path.abspath(filename)[len(os.path.abspath(env.shared_resource)):-3])
                if os.path.isfile(new_filename) and not force:
                    if md5 is not None and md5 != calculateMD5(new_filename, partial=True):
                        env.error('MD5 signature mismatch: {} (signature: {} calculated: {})'
                            .format(new_filename, md5, calculateMD5(new_filename, partial=True)))
                    else:
                        env.log('Reusing existing decompressed file {}'.format(new_filename))
                        return new_filename
            else:
                raise RuntimeError('Failed to decompress file {}: directory not writable'.format(filename))
        # dest_dir can be '' if there is no path for filename
        if dest_dir and not os.path.isdir(dest_dir):
            os.makedirs(dest_dir)
        #
        env.log('Decompressing {} to {}'.format(filename, new_filename))
        try:
            with gzip.open(filename, 'rb') as input, open(new_filename, 'wb') as output:
                buffer = input.read(100000)
                while buffer:
                    output.write(buffer)
                    buffer = input.read(100000)
            if md5 is not None and md5 != calculateMD5(new_filename, partial=True):
                env.error('MD5 signature mismatch: {} (signature {}, calculated {})'
                    .format(new_filename, md5, calculateMD5(new_filename, partial=True)))
        # Python 2.7.4 and 3.3.1 have a regression bug that prevents us from opening
        # certain types of gzip file (http://bugs.python.org/issue17666).
        except TypeError as e:
            raise RuntimeError('Failed to open gzipped file {} due to a bug '
                'in Python 2.7.4 and 3.3.1. Please use a different '
                'version of Python or decompress this file manually.'.format(filename))
        #
        if inplace:
            try:
                os.remove(filename)
            except:
                pass
        return new_filename
    else:
        return filename

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

def physicalMemory():
    '''Get the amount of physical memory in the system'''
    # MacOSX?
    if platform.platform().startswith('Darwin'):
        # FIXME
        return None
    elif platform.platform().startswith('Linux'):
        try:
            res = subprocess.check_output('free').decode().split('\n')
            return int(res[1].split()[1])
        except Exception as e:
            return None

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

#http://stackoverflow.com/a/13197763, by Brian M. Hunt
class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = newPath

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

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
            if sys.version_info.major == 3:
                instream = instream.encode(sys.getdefaultencoding())
            out, error = tc.communicate(instream)
        else:
            out, error = tc.communicate()
        if sys.version_info.major == 3:
            out = out.decode(sys.getdefaultencoding())
            error = error.decode(sys.getdefaultencoding())
        if return_zero:
            if tc.returncode < 0:
                raise ValueError ("Command '{0}' was terminated by signal {1}".format(cmd, -tc.returncode))
            elif tc.returncode > 0:
                raise ValueError ("{0}".format(error))
        if error.strip() and show_stderr:
            env.log(error)
    except OSError as e:
        raise OSError ("Execution of command '{0}' failed: {1}".format(cmd, e))
    # everything is OK
    if upon_succ:
        # call the function (upon_succ) using others as parameters.
        upon_succ[0](*(upon_succ[1:]))
    return out.strip(), error.strip()

def zipdir(path, zipfile, arcroot = '/'):
    path = os.path.normpath(path)
    for root, dirs, files in os.walk(path):
        for f in files:
            zipfile.write(os.path.join(root, f), arcname = os.path.join(arcroot, root[len(path) + 1:], f))

def removeFiles(dest, exclude = [], hidden = False):
    if os.path.isdir(dest):
        for item in os.listdir(dest):
            if item.startswith('.') and hidden == False:
                continue
            if os.path.splitext(item)[1] not in exclude:
                try:
                    os.remove(os.path.join(dest,item))
                except:
                    pass

def removeEmptyDir(directory):
    try:
        if not os.listdir(directory):
            os.rmdir(directory)
    except:
        pass

def copyFiles(pattern, dist, ignore_hidden = True):
    mkpath(dist)
    for fl in glob.glob(pattern):
        if os.path.isfile(fl):
            shutil.copy(fl, dist)

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
        if isinstance(v, MutableMapping):
            items.extend(flatten_dict(v).items())
        else:
            items.append((k, v))
    return dict(items)
