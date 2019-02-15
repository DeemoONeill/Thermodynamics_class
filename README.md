# Thermodynamics_class

Contains a class for calculating DG, Keq, and ultimately reactant and product concentrations for 
thermodynamically limited reactions (DG > 0).

Example use cases:


Single calculation:
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
            
Multiple temperatures:
----------------------
        >>> DMC = Thermo(DH=-24, DS=-0.123)
        >>> DMC.Reactants(["CO2", "MeOH"], [8, 24], [1, 2])
        >>> DMC.Products(["DMC", "H2O"], [0, 0], [1, 1])
        >>> T = 298
        >>> for i in range(2)
        ...     DMC.thermodynamics(T)
        ...     T += 5
        >>> print(DMC.return_dictionary())
        {'DH': [-24, -24], 'DS': [-0.123, -0.123], 'T': [298, 303],
        'CO2': [4.921006688684159, 5.067904590645182],
        'MeOH': [17.842013377368318, 18.13580918129034],
        'DMC': [3.0789933113158408, 2.9320954093548233],
        'H2O': [3.0789933113158408, 2.9320954093548233]}

Multiple starting concentrations:
---------------------------------
        >>> DMC = Thermo(DH=-24, DS=-0.123)
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

Product starting concentrations:
--------------------------------
        >>> DMC = Thermo(DH=-24, DS=-0.123)
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
