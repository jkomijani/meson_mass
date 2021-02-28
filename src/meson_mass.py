#!/usr/bin/env python3

# Created by Javad Komijani, University of Tehran, on 21/Feb/2021.
# Copyright (C) 2021 Javad Komijani

"""
To describe the mass of a generic heavy-light pseudoscalar meson
simulated on a lattice with lattice spacing *a* and volume *V*,
the Fermilab Lattice, MILC, and TUMQCD uses a function as

.. math::

    M_{H_x} (\{m_h, m_x\}, \{m'_l, m'_l, m'_s, m'_c\}, \{a, V\};
             \{p_1,\cdots,p_{67}\}),

which takes as inputs three sets of values and contains 67 parameters.
The inputs are the simulation masses of the valence heavy and light quarks
(in the first braces), the simulation masses of four flavors of sea quarks
(in the second braces), and the lattice spacing and volume (in the third
braces).
Note that the lattice simulations are usually performed at the isospin limit,
where the masses of two light quarks (up and down quarks) are equal.
We use the term *light* quark to call these two quarks at the isospin limit
and :math:`m_l` to denote their masses.

The 67 parameters of :math:`M_{H_x}` are determined statistically by fitting
the function to lattice-QCD data.
The continuum and infinite-volume limits of the function are then obtained
by setting :math:`a = 0` and :math:`V = \infty`.

One can use the function :math:`M_{H_x}` to explore how a heavy-light meson
mass depends on the masses of its components.
For applications in ChPT, one may be interested to fix the heavy and strange
quark masses and vary the light quark mass.
This reduces the function :math:`M_{H_x}` to a one-variable function,
which is still relatively complicated because of non-analytic terms stemming
from chiral perturbation theory; see Ref1_ and Ref2_ for details.
(Note that the non-analytic terms are very small in practice.)
When quark masses are neither too big nor too small, the nonanalytic terms
can be replaced with polynomials.

This module provides four (cubic) polynomial functions to describe the masses
of :math:`D, D_s, B` and :math:`B_s` mesons in terms of the light quark mass
in units of the strange quark mass, i.e. in terms of :math:`x = m_l/m_s`.
(As said above, the strange, charm, and bottom masses are set to their physical
values.)
The functions take :math:`x = m_l/m_s` as input and return the mass of
the corresponding mesons (in MeV) and their one-sigma uncertainties.
The pacakge `gvar`_ is used to handle the values and their errors.
If an array of light quark masses (in units of strange mass) are provided,
each function returns an array of `gvar`_ variables that contains
the correlation between variables.

.. _gvar: https://github.com/gplepage/gvar

The following figure illustrates the :math:`D, D_s, B` and :math:`B_s` meson
masses against :math:`x = m_l/m_s`:

.. image:: observation_vs_polyfit.jpeg

The error bars illustrate the results obtained from the analysis of the project
on Ref1_, and the error bands are derived from the effective formulas of this
package. The domain of validity of the functions are restricted to
:math:`x\in[0.02, 0.5]`.

For the :math:`D` system, e.g. with quark mass ratios 0.1, 0.2, and
0.3, one can simply use:

    >>> from meson_mass import Model
    >>> Model('D').predict([0.1, 0.2, 0.3])
    >>> array([1876.81(97), 1890.8(1.3), 1904.6(1.5)], dtype=object)

where output is in MeV and provided as an array of gvar_ variables.
(One can use `gvar.mean(...)` `gvar.sdev(...)` to access the mean and one-sigma
error in the output.)
Similarly, for :math:`D_s, B` and :math:`B_s` systems, we have

    >>> Model('Ds').predict([0.1, 0.2, 0.3])
    >>> array([1968.87(25), 1971.87(46), 1975.11(94)], dtype=object)

    >>> Model('B').predict([0.1, 0.2, 0.3])
    >>> array([5288.2(1.3), 5302.7(1.6), 5316.9(2.0)], dtype=object)

    >>> Model('Bs').predict([0.1, 0.2, 0.3])
    >>> array([5369.82(41), 5374.25(83), 5378.9(1.4)], dtype=object)

Command line inteface is available too::

    $ meson_mass.py --model D  0.1 0.2 0.3
    $ [1876.81(97) 1890.8(1.3) 1904.6(1.5)]
    $
    $ meson_mass.py --model Ds  0.1 0.2 0.3
    $ [1968.87(25) 1971.87(46) 1975.11(94)]
    $
    $ meson_mass.py --model B  0.1 0.2 0.3
    $ [5288.2(1.3) 5302.7(1.6) 5316.9(2.0)]
    $
    $ meson_mass.py --model Bs  0.1 0.2 0.3
    $ [5369.82(41) 5374.25(83) 5378.9(1.4)]


The parameters of this module correspond to a QCD-only world;
in lattice simulations of heavy-light mesons, the QED effects are not included.
To tune the charm and bottom quark masses, the analysis uses the experimental
values of :math:`D_s` and :math:`B_s` masses, but it subtracts their
QED contributions in a specific scheme.
To compare the results of this module with experimental values, which include
QED effects, an option is provided to adjust the outputs depending on
the charges of the heavy and light components of the meson.
For the masses of :math:`D_s` and :math:`B_s` mesons, we have

    >>> m_ud = 0.0368  # average of u and d mases in unit of the strange mass
    >>> Model('Ds', QCD_only=False, e_heavy=2/3, e_light=1/3).predict(ml)
    >>> 1968.27(10)
    >>> Model('Bs', QCD_only=False, e_heavy=-1/3, e_light=1/3).predict(ml)
    >>> 5366.82(22)

These results are equal to the experimental masses
:math:`M_{D_s} = 1968.27(10)` and :math:`M_{B_s} = 5366.82(22)`
that are used as inputs to fix the charm and bottom quark masses, respectively.
(It is not surprising that this module gives such a precise value for these two
quantities.)

As a concluding remark, it should be emphasized that the effective formula
are expected to be used at :math:`x\in[0.02, 0.5]`.
The results of the code outside the given domain cannot be trusted.
For instance, evaluating the :math:`D` mass at the chiral limit, one obtains::

    $ meson_mass.py --model D 0
    $ 1862.45(38)

while the *correct* value reported in eq (5.31) Ref1_ is 1862.3(1.3) MeV.
This rapid increase of uncertainty indicates the importance of non-analytic
terms in approaching the chiral limit.

.. _Ref1: https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.054517
.. _Ref2: https://journals.aps.org/prd/abstract/10.1103/PhysRevD.97.034503

"""

import numpy as np
import gvar as gv


class Model:
    """..."""

    def __init__(self, meson_name=None, **model_kwargs):

        self.param = Param(meson_name, **model_kwargs)

    def fit(self):
        """ An abstract method."""

    def predict(self, x):
        """ Return the mass of a heavy-light meson in MeV.

        Note that the effective formula are originally organized as functions
        of :math:`z = m_l/ (0.4 m_s)`, but the interface of function uses
        :math:`x = m_l/m_s` as input.

        Parameters:
            x : array-like or scalar
                Ratio of the light quark mass to the strange quark mass,
                :math:`x = m_l/m_s`.

            par : an instance of class `par`
                The class of coefficients of polynomial expansion in
                powers of (0.4 times) :math:`x`.

        Returns:
            mass : gvar
                The mass(es) of the meson at given `x` as a gvar object.
                To get the mean value of the meson mass(es) use `gv.gvar(mass)`
                and to get its one-sigma uncertainty use `gv.sdev(mass)`.
        """

        if np.ndim(x) == 1:
            x_ = x
        elif np.ndim(x) == 0:
            x_ = [x]
        else:
            raise Exception("x is expected to be a scalar or 1-dim array like")

        z = 2.5 * np.array(x_)  # converting from `s` units to `p4s` units

        par = self.param
        Z = np.vander(z, par.degree+1, increasing=True)
        mean = np.sum(Z * par.mean, axis=1)
        cov = Z @ par.cov @ Z.transpose()

        if np.ndim(x) == 0:
            return gv.gvar(mean[0], cov[0, 0]**0.5)  # gv.gvar(mean, sigma)
        else:
            return gv.gvar(mean, cov)  # gv.gvar(mean, covariance_matrix)

    def set_param(self, **kwargs):
        """..."""

        self.param.__dict__.update(**kwargs)


class Param:
    """ This is a class for parameters of effective formula describing
    the masses of heavy-light mesons.

    Four systems are considered in this class: :math:`D, D_s, B` and
    :math:`B_s` systems.

    Attributes:
        degree : int
            The degree of the polynomial used for the effective formula.

        mean : array
            The mean value of the coefficients of the polynomial expansion.
            The coefficients are in MeV and the expansion parameter is
            dimensionless.

        sigma : array
            One sigma uncertainty of the coefficients of the expansion.

        corr : 2-dim array
            The correlation matrix of the coefficients of the expansion.

        cov : 2-dim array
            The covariance matrix of the coefficients of the expansion.

    """

    def __init__(self, meson_name=None, QCD_only=True, **electric_charges):
        """..."""
        if meson_name in ('D', 'charm-light'):
            self._charm_light()
        elif meson_name in ('Ds', 'charm-strange'):
            self._charm_strange()
        elif meson_name in ('B', 'bottom-light'):
            self._bottom_light()
        elif meson_name in ('Bs', 'bottom-strange'):
            self._bottom_strange()

        if not QCD_only:
            self.qed_shift(**electric_charges)

    def qed_shift(self, e_heavy=None, e_light=None):
        """Adjust the meson masses for QED contributions.

        The parameters of this module correspond to a QCD-only world,
        with a specific scheme to convert the real world results to the
        QCD-only world.

        This method allows to adjust the parameters for QED contributions.
        To this end, we use the relations of the `paper`_.

        Parameters
        ----------
        e_heavy : float
            The charge of the heavy component (quark/antiquark)
        e_light : float
            The charge of the light component (quark/antiquark)
        """

        A = 4.44 # MeV  (EM contribution to the potential from a Coulomb type interaction)
        B = 2.43 # MeV  (EM contribution to the self-energy of the light quark)
        
        qed_shift = A*e_heavy*e_light + B*e_light**2

        self.mean[0] += qed_shift

    def _charm_light(self):
        self.degree = 3
        self.mean = np.array([1862.45099286, 58.68981998, -5.74628366, 3.20775425])  # MeV
        self.sigma = np.array([0.38477103, 4.63868325, 6.33864591, 2.38286455])
        self.corr = np.array(
            [[ 1.        ,  0.22559089, -0.33935846,  0.51497625],
             [ 0.22559089,  1.        , -0.94273073,  0.7063392 ],
             [-0.33935846, -0.94273073,  1.        , -0.74129706],
             [ 0.51497625,  0.7063392 , -0.74129706,  1.        ]]
            )
        # corr_eigenvalues = np.array([2.80961397, 0.87849176, 0.26213027, 0.04976401])
        self.cov = np.diag(self.sigma) @ self.corr @ np.diag(self.sigma)

    def _charm_strange(self):
        self.degree = 3
        self.mean = np.array([1965.93050, 11.8318531, -0.806120707, 1.80387062])  # MeV
        self.sigma = np.array([0.22105652, 2.42501362, 3.61272179, 1.29586647])
        self.corr  = np.array(
            [[ 1.        , -0.88971983,  0.74939704,  0.08170989],
             [-0.88971983,  1.        , -0.87476987, -0.0696281 ],
             [ 0.74939704, -0.87476987,  1.        , -0.127539  ],
             [ 0.08170989, -0.0696281 , -0.127539  ,  1.        ]]
            )
        # corr_eigenvalues = np.array([2.67759947, 1.03391197, 0.22265057, 0.06583799])
        self.cov = np.diag(self.sigma) @ self.corr @ np.diag(self.sigma)

    def _bottom_light(self):
        self.degree = 3
        self.mean = np.array([5273.06246, 62.2177897, -7.67618280, 3.64580249])  # MeV
        self.sigma = np.array([1.38390992, 6.66712555, 7.77065485, 2.89792283])
        self.corr  = np.array(
            [[ 1.        , -0.3960323 ,  0.00183392,  0.05570996],
             [-0.3960323 ,  1.        , -0.87568809,  0.72434404],
             [ 0.00183392, -0.87568809,  1.        , -0.8290985 ],
             [ 0.05570996,  0.72434404, -0.8290985 ,  1.        ]])
        # corr_eigenvalues = np.array([2.64532559, 1.11231589, 0.2028819 , 0.03947662])
        self.cov = np.diag(self.sigma) @ self.corr @ np.diag(self.sigma)

    def _bottom_strange(self):
        self.degree = 3
        self.mean = np.array([5365.41675, 17.7713126, -1.11310610, 1.81987915])  # MeV
        self.sigma = np.array([0.3280838 , 2.86605657, 3.61171569, 1.30380341])
        self.corr  = np.array(
            [[ 1.        , -0.73926191,  0.48106424,  0.08698118],
             [-0.73926191,  1.        , -0.7081466 , -0.0980904 ],
             [ 0.48106424, -0.7081466 ,  1.        , -0.13903376],
             [ 0.08698118, -0.0980904 , -0.13903376,  1.        ]]
            )
        # corr_eigenvalues = np.array([2.2928649 , 1.05889039, 0.47118009, 0.17706462])
        self.cov = np.diag(self.sigma) @ self.corr @ np.diag(self.sigma)


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    MESON_METAVAR = "{D, Ds, B, Bs, charm-light, ...}"
    parser.add_argument('inputs', metavar='input', type=float, nargs='+')
    parser.add_argument("--model",  dest="meson_name", metavar=MESON_METAVAR,
                        type=str, default=None)
    try:
        args = parser.parse_args()
        x_list = args.inputs if len(args.inputs)>1 else float(args.inputs[0])
        print(Model(args.meson_name).predict(x_list))
    except (IndexError, ValueError, AttributeError):
        print("usage, e.g.: meson_mass.py --model D  0.1  0.15 ...")
        print("supported models: ", MESON_METAVAR)
