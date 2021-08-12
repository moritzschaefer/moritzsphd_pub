import logging
import os
import re
from sqlite3 import OperationalError

import gffutils
import numpy as np
import pandas as pd
from sklearn.compose import \
    TransformedTargetRegressor as sklearnTransformedTargetRegressor

from moritzsphd.util.file import DATA_DIR, hash, remote_file


def average_sample(df):
    '''
    Return the average of replicates
    '''
    return pd.concat([
        df.select_dtypes(exclude=[np.number]), df.select_dtypes(include=[np.number]).groupby(
            by=lambda f: f[:f.rfind('_')], axis=1).mean()
    ],
                     axis=1)


def generate_gffutils_db(gff_filename):
    db_filename = DATA_DIR / hash(str(gff_filename))
    try:
        return gffutils.create_db(
            str(gff_filename), str(db_filename),
            merge_strategy='merge')  # might be an issue, alternative:  'warning'
    except OperationalError as e:
        try:
            return gffutils.FeatureDB(db_filename)
        except TypeError:  # corrupted file
            os.remove(db_filename)
            raise RuntimeError('DB or GFF Files corrupted. Deleted db')
            #os.remove(gff_filename)
            #return gffutils.create_db(
            #    str(gff_filename), str(db_filename),
            #    merge_strategy='merge')  # might be an issue, alternative:  'warning'


def parse_attributes(gene, index=-1):
    '''
    parse attributes as found in a gtf file
    index
    '''
    try:
        if '=' in gene[index]:
            return dict(v.split('=') for v in gene[index].split(';') if '=' in v)
        else:
            return dict(re.match('([^ ]*) "(.*)"', v).groups() for v in gene[index].strip(';').split('; ') if ' ' in v)
    except KeyError:
        return {}


def log_call(func, log_level=logging.INFO):
    def wrapper(*args, **kwargs):
        logging.log(log_level, f'Entering {func.__name__}')
        func(*args, **kwargs)
        logging.log(log_level, f'Leaving {func.__name__}')

    return wrapper


# https://stackoverflow.com/questions/58155778/target-transformation-and-feature-selection-in-scikit-learn
class TransformedTargetRegressor(sklearnTransformedTargetRegressor):
    @property
    def feature_importances_(self):
        try:
            return self.regressor_.feature_importances_
        except AttributeError:  # maybe we are in a pipeline
            return self.regressor_.steps[-1][-1].feature_importances_

    @property
    def coef_(self):
        try:
            return self.regressor_.coef_
        except AttributeError:  # maybe we are in a pipeline
            return self.regressor_.steps[-1][-1].coef_
