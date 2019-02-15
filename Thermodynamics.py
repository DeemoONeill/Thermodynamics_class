# -*- coding: utf-8 -*-
"""
Created on Tue Feb  5 16:22:54 2019

@author: Matt O'Neill
"""

import math as ma


class Thermo:
    """Takes values of DH and DS at a given temperature
    ------------------
    Keyword arguments:
    ------------------
    DH -- Enthalpy change of the reaction at Temperature T (float)
    DS -- Entropy change of the reaction at Temperature T (float)
    T -- Temperature at which DH and DS were calculated (default 298)

    ------------
    Key Methods:
    ------------
    Reactants(names, starting_concentrations, stoichiometry)
    Products(names, starting_concentrations, stoichiometry)
    sets properties of the reactants and products - names,
    starting_concentrations and stoichiometrys should all be given as lists
    (even for single values)

    dictionary_flush()
    Sets all keys in the dictionary to empty lists

    return_dictionary()
    Returns the dictionary of calculated values includes T, DH, DS, product
    and reactant concentrations

    thermodynamics(temperature = 298, error=1e-09, p=20, K=True,
    return_dict = False)
    calculates reactant and product concentrations for the given temperature.
    Takes Temperature, error, P, K, and return_dict arguments.
    error gives the tolerance between Keq calculated from DG and Keq calculated
    by products/reactants
    p is used in the same was as in PID systems, it is multiplied by the
    difference between the setpoint and current measurement, and used to adjust
    the concentration values.
    K is whether the values for temperature are saved in Kelvin or Celsius
    (no support for F)
    return_dict returns the dictionary at the end of the calculation. mainly
    used for single     calculations
    """
    def __init__(self, DH=None, DS=None, T=298):
        self.DSstd = DS
        self.DHstd = DH
        self._DS_T = self.DSstd
        self._DH_T = self.DHstd
        self._dictionary = {'DH': [], 'DS': [], 'T': []}
        self._T0 = T

    def _Calc_DG(self, T):
        self.DG = self._DH_T - T * self._DS_T
        return self.DG

    def Reactants(self, names, starting_concentrations, stoichiometry):
        """Sets the reactant properties for the reaction
        ------------------
        Keyword arguments:
        ------------------
        names -- a list of chemical names
        starting_concentrations -- a list of starting concentrations
        stoichiometry -- takes a list of the reaction stoichiometry
        --------
        example:
        --------
        for the reaction CO2 + 2MeOH --> DMC + H2O

        dmc = thermo(DH=-24, DS=-0.12, T=298)
        dmc.Reactants(['CO2', 'MeOH'], [8, 24.5], [1,2])
        """
        self._ReactantNames = names
        self._ReactantStartingConcentration = starting_concentrations
        self._ReactantStoichiometry = stoichiometry
        # self._ReactantHeatCapacities = heat_capacities
        for i in range(len(self._ReactantNames)):
            try:
                self._dictionary[self._ReactantNames[i]]
            except KeyError:
                self._dictionary[self._ReactantNames[i]] = []

    def Products(self, names, starting_concentrations, stoichiometry,
                 heat_capacities=None):
        """Sets the product properties for the reaction
        ------------------
        Keyword arguments:
        ------------------
        names -- a list of chemical names
        starting_concentrations -- a list of starting concentrations
        stoichiometry -- takes a list of the reaction stoichiometry
        --------
        example:
        --------
        for the reaction CO2 + 2MeOH --> DMC + H2O

        dmc = thermo(DH = -24, DS = -0.12, T = 298)
        dmc.Products(['DMC', 'H2O'], [0, 0], [1,1])

        """
        self._ProductNames = names
        self._ProductStartingConcentration = starting_concentrations
        self._ProductStoichiometry = stoichiometry
        self._ProductHeatCapacities = heat_capacities
        for i in range(len(self._ProductNames)):
            try:
                self._dictionary[self._ProductNames[i]]
            except KeyError:
                self._dictionary[self._ProductNames[i]] = []

    def _nCP(self):
        heat_capacity = 0
        for i in range(len(self._ProductStoichiometry)):
            heat_capacity += (self._ProductStoichiometry[i] *
                              self._ProductHeatCapacities[i])
        for i in range(len(self._ReactantStartingConcentration)):
            heat_capacity -= (self._ReactantStoichiometry[i] *
                              self._ReactantHeatCapacities[i])
        return heat_capacity

    def _correction(self, T):
        self._DH_T = self.DHstd - self._nCP()*(T-self._T0)
        self._DS_T = self.DSstd - self._nCP()*ma.log(T/self._T0)

    def _appender(self, T):
        self._dictionary['DS'].append(self._DS_T)
        self._dictionary['DH'].append(self._DH_T)
        self._dictionary['T'].append(T)

    def flush_dictionary(self):
        """Flushes the dictionary, sets each key to an empty list
        """
        for key in self._dictionary:
            self._dictionary[key] = []

    def return_dictionary(self):
        """Returns the data containing dictionary"""
        return self._dictionary

    def thermodynamics(self, temperature=298, error=1e-09, p=20, K=True,
                       return_dict=False):
        """Calculates the equilibrium concentrations for each of the products
        and reactants in the reaction at a given temperature.

        first DG is calculated using:
        DG = DH - TDS

        followed by the equilibrium constant Keq:
        Keq = exp(-DG/RT)

        The concentrations of each of the reactants and products are then
        adjusted until this valuse of Keq is reached using the equation:
        Keq = [Products]**S/[Reactants]**S

        ------------------
        Keyword arguments:
        ------------------

        temperature -- The temperature in Kelvin at which the equilibrium
        concentrations should be calculated (default 298)

        error -- tolerance for when Keq calculated from DG and from
        concentrations are considered equal (default 1e-9)

        p -- The proportional correction applied to the result. if calculation
        is slow, try adjusting the value for P (default 20)

        K -- whether or not Tempertature is appended to the dictionary in
        Kelvin (True) or °C (False). (default True)

        return_dict -- whether to return the dictionary after calculation
        (default False)

        ---------
        examples:
        ---------
        -------------------
        single calculation:
        -------------------
        >>> DMC = Thermo(DH=-24, DS=-0.123)
        >>> DMC.Reactants(["CO2", "MeOH"], [8, 24], [1, 2])
        >>> DMC.Products(["DMC", "H2O"], [0, 0], [1, 1])
        >>> DMC.thermodynamics(298)
        >>> print(DMC.return_dictionary())
            {'DH': [-24], 'DS': [-0.123], 'T': [298],
            'CO2': [4.921006688684159],
            'MeOH': [17.842013377368318],
            'DMC': [3.0789933113158408],
            'H2O': [3.0789933113158408]}

        ----------------------
        multiple temperatures:
        ----------------------

        >>> DMC = Thermo(DH= -24, DS= -0.123)
        >>> DMC.Reactants(["CO2", "MeOH"], [8, 24], [1, 2])
        >>> DMC.Products(["DMC", "H2O"], [0, 0], [1, 1])
        >>> T = 298
        >>> for i in range(2):
        ...     DMC.thermodynamics(T)
        ...     T += 5

        >>> print(DMC.return_dictionary())
        {'DH': [-24, -24], 'DS': [-0.123, -0.123], 'T': [298, 303],
        'CO2': [4.921006688684159, 5.067904590645182],
        'MeOH': [17.842013377368318, 18.13580918129034],
        'DMC': [3.0789933113158408, 2.9320954093548233],
        'H2O': [3.0789933113158408, 2.9320954093548233]}

        ---------------------------------
        Multiple starting concentrations:
        ---------------------------------

        >>> DMC = Thermo(DH= -24, DS= -0.123)
        >>> DMC.Reactants(["CO2", "MeOH"], [8, 24], [1, 2])
        >>> DMC.Products(["DMC", "H2O"], [0, 0], [1, 1])
        >>> T = 298
        >>> CO2 = 8
        >>> MeOH = 10
        >>> for i in range(2):
        ...     DMC.Reactants(["CO2", "MeOH"], [CO2, MeOH], [1, 2])
        ...     DMC.thermodynamics(T, return_dict = True)
        ...     CO2+=5
        >>> print(DMC.return_dictionary())
        {'DH': [-24, -24], 'DS': [-0.123, -0.123], 'T': [298, 298],
        'CO2': [6.574179093358968, 11.283789139423279],
        'MeOH': [7.148358186717942, 6.56757827884657],
        'DMC': [1.4258209066410286, 1.7162108605767157],
        'H2O': [1.4258209066410286, 1.7162108605767157]}

        --------------------------------
        Product starting concentrations:
        --------------------------------

        >>> DMC = Thermo(DH= -24, DS= -0.123)
        >>> DMC.Reactants(["CO2", "MeOH"], [8, 24], [1, 2])
        >>> DMC.Products(["DMC", "H2O"], [0, 0.134], [1, 1])
        >>> T = 298
        >>> CO2 = 8
        >>> MeOH = 10
        >>> for i in range(2):
        ...     DMC.Reactants(["CO2", "MeOH"], [CO2, MeOH], [1, 2])
        ...     DMC.thermodynamics(T, return_dict = True)
        ...     CO2+=5
        >>> print(DMC.return_dictionary())
        {'DH': [-24, -24], 'DS': [-0.123, -0.123], 'T': [298, 298],
        'CO2': [6.6175669997054625, 11.324869902582156],
        'MeOH': [7.235133999410933, 6.649739805164333],
        'DMC': [1.3824330002945329, 1.6751300974178351],
        'H2O': [1.5164330002945328, 1.809130097417835]}

        """

        self._Calc_DG(temperature)
        DGkeq = ma.exp(-(self.DG*1000)/(8.314*temperature))
        Conckeq = 0
        pconc = self._ProductStartingConcentration.copy()
        sconc = self._ReactantStartingConcentration.copy()

        while ma.isclose(DGkeq, Conckeq, rel_tol=error) is False:
            iterator = p*(DGkeq-Conckeq)
            numerator = 1
            denominator = 1
            # Reactants
            for i in range(len(sconc)):
                working = sconc.pop(0)
                working -= iterator*self._ReactantStoichiometry[i]
                denominator *= working ** self._ReactantStoichiometry[i]
                sconc.append(working)

            # Products
            for i in range(len(pconc)):
                working = pconc.pop(0)
                working += iterator*self._ProductStoichiometry[i]
                numerator *= working ** self._ProductStoichiometry[i]
                pconc.append(working)

            Conckeq = numerator / denominator

        for i in range(len(sconc)):
            self._dictionary[self._ReactantNames[i]].append(sconc[i])
        for i in range(len(pconc)):
            self._dictionary[self._ProductNames[i]].append(pconc[i])
        if K:
            self._appender(temperature)
        else:
            self._appender((temperature-273))
        if return_dict:
            self.return_dictionary()
