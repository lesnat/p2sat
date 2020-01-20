#coding:utf8
import numpy as np

class _Dataset(object):
    r"""
    Mother class for all fp2sat datasets.
    """
    def __len__(self):
        r"""
        Return the total number of lines of the dataset.
        """
        return len(self.read.w)
