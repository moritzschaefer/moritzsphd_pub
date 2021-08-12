import gzip
import hashlib
import logging
import os
import shutil
import socket
import subprocess
import tarfile
import tempfile
import zipfile
from ftplib import FTP
from pathlib import Path
from urllib.parse import urlparse

import requests

try:
    from xdg.BaseDirectory import save_cache_path
    DATA_DIR = Path(save_cache_path('moritzsphd'))
    if not DATA_DIR.exists():
        os.mkdir(DATA_DIR)
except ImportError:
    DATA_DIR = Path('/home/moritz/.cache/moritzsphd/')


def hash(s):
    try:
        return hashlib.sha256(s).hexdigest()
    except TypeError:
        return hashlib.sha256(s.encode('utf-8')).hexdigest()


def delete_cached_file(url):
    '''
    Helper function to delete an erroneous file or something like that..
    '''
    filename = DATA_DIR / hash(url)
    try:
        os.remove(filename)
    except FileNotFoundError:
        logging.warn(f'File with URL {url} was not fount locally ({filename})')


def ssh_file(url):
    tf = tempfile.NamedTemporaryFile(delete=False)
    tf.close()
    bash_args = ['scp', str(url), str(tf.name)]
    p = subprocess.run(bash_args)

    if p.returncode != 0:
        raise RuntimeError(' '.join(bash_args), p.returncode, p.stderr, p.stdout)

    return tf


def http_file(url):
    r = requests.get(url, allow_redirects=True)

    tf = tempfile.NamedTemporaryFile(delete=False)
    tf.write(r.content)
    tf.close()
    return tf


def ftp_file(url):
    tf = tempfile.NamedTemporaryFile(delete=False)
    parsed_url = urlparse(url)
    ftp = FTP(parsed_url.netloc)
    ftp.login('anonymous', '')

    ftp.retrbinary("RETR " + parsed_url.path, tf.write)
    tf.close()
    return tf

# TODO use fsspec python library for caching
def remote_file(url, zip_extract_name=None, gunzip=False, tar_extract_name=None, reload=False, with_extension=False):
    '''
    Automatically downloads a URL (if not cached) and provides the file path;
    Also tries to automatically unzip files (only works with ZIP files containing a single file with the correct naming..)
    :zip_extract_name: extract the specified filename from a ZIP file
    :gunzip: gzip extract the filename after downloading
    :tar_extract_name: same as for zip_extract_name (but for tar files)
    '''

    if url.startswith('cclab:') and socket.gethostname() == 'mhs-cclab-srv001':
        return Path(url.lstrip('cclab:'))

    if with_extension:
        filename = DATA_DIR / (hash(url) + os.path.splitext(url)[1])
    else:
        filename = DATA_DIR / hash(url)

    if filename.exists() and reload:
        os.remove(filename)

    if not filename.exists():
        if url.startswith('http'):
            tf = http_file(url)
        elif url.startswith('ftp'):
            tf = ftp_file(url)
        else:
            tf = ssh_file(url)


        with open(filename, 'wb') as f_out:
            if tar_extract_name:
                with tarfile.open(tf.name, 'r:' if not gunzip else 'r:gz') as tared_file:
                    f_out.write(tared_file.extractfile(tar_extract_name).read())
            elif zip_extract_name:
                with zipfile.ZipFile(tf.name, 'r') as zipped_file:
                    f_out.write(zipped_file.read(zip_extract_name))
            else:
                if gunzip:
                    f_in = gzip.open(tf.name, 'rb')
                else:
                    f_in = open(tf.name, 'rb')

                shutil.copyfileobj(f_in, f_out)
                f_in.close()
        os.remove(tf.name)

    return filename
