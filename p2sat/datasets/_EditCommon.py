#coding:utf8
import numpy as np

class _EditCommon(object):
    r"""
    """
    def update(self):
        raise NotImplementedError()

    def generate(self):
        raise NotImplementedError()

    def normalize(self, norm=1., verbose=True):
        r"""
        Normalize the weights to a given value.

        Parameters
        ----------
        norm : float, optional
            Total weight. Default is 1.
        verbose : bool, optional
            verbosity
        """
        if verbose:print("Normalizing weights to %.2E ..."%norm)

        w = norm * self._ds.read.w/sum(self._ds.read.w)

        if verbose:print("Done !")

        self.update(w, *self._ds.read.dataset[1:], in_code_units=False, verbose=verbose)

    def round_quantity(self, qty, decimals=8,verbose=True):
        r"""
        Round the given quantity.

        Parameters
        ----------
        qty : str
            quantity to round
        decimals : int, optional
            number of decimals
        verbose : bool, optional
            verbosity

        Examples
        --------
        >>> eps = ExamplePhaseSpace()
        >>> eps.read.x
        array([-4.26043957, -8.225104  ,  0.25424565, ..., -3.19180518])

        >>> eps.edit.round_quantity("x", decimals=1)
        Rounding qty x ...
        Done !

        >>> eps.read.x
        array([-4.3, -8.2,  0.3, ..., -3.2])
        """
        if verbose:print("Rounding quantity %s ..."%qty)
        qty = self._ds.read.quantity(qty)
        dataset = []
        for q in self._ds.read.dataset:
            if q is qty:
                dataset.append(np.around(qty, decimals))
            else:
                dataset.append(q)

        if verbose: print("Done !")
        self.update(*dataset, verbose=verbose)
