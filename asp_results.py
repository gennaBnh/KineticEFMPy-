import pickle


#-------------CODE FROM THE "read_output.ipynb" FILE OF THE ASP-EFM LIBRARY--------------

def recup_efms(efmFile):
    
    with open(efmFile, 'rb') as f:
        efms = pickle.load(f)

    # select only biomass-producing efms
    df_biomass = efms.copy()
    df_biomass = df_biomass.loc[df_biomass['biomass'] > 0]

    # compute operating costs or yields
    # units are O2 Mol / C Mol Biomass and C Mol Glucose / C Mol Biomass
    # 1 Glucose molecule corresponds to 6 C Mol
    df_yields = df_biomass.copy()
    df_yields['o2_biomass_yield'] = df_yields['ex-o2'] / df_yields['biomass']
    df_yields['c_biomass_yield'] = 6*df_yields['pts'] / df_yields['biomass']
    df_yields
    num_c_moles_biomass = 6*106.86
    df_yields['o2_biomass_yield'] /= num_c_moles_biomass
    df_yields['c_biomass_yield'] /= num_c_moles_biomass
    df_yields

    # filter efms with yields that are outliers
    df_yields_filter = df_yields.copy()
    df_yields_filter = df_yields_filter.loc[df_yields['o2_biomass_yield'] < 1]
    df_yields_filter = df_yields_filter.loc[df_yields['c_biomass_yield'] < 10]
    df_yields_filter

    # returns True if the row is at the pareto frontier for variables xlabel and ylabel
    config_base = df_yields_filter
    def is_pareto_front(row, xlabel, ylabel):
        
        x = row[xlabel]
        y = row[ylabel]
        
        # look for points with the same y value but smaller x value
        is_min_x = config_base.loc[config_base[ylabel]==y].min()[xlabel] >= x
        # look for points with the same x value but smaller y value
        is_min_y = config_base.loc[config_base[xlabel]==x].min()[ylabel] >= y
        # look for points that are smaller in both x and y
        is_double = len(config_base.loc[(config_base[xlabel]<x) & (config_base[ylabel]<y)])==0
        
        return is_min_x and is_min_y and is_double

    # array of True/False indicating whether the corresponding row is on the pareto frontier
    is_pareto = df_yields.apply(lambda row: is_pareto_front(row, 'c_biomass_yield', 'o2_biomass_yield'), axis=1)

    # compute pareto optimal efms
    is_pareto[is_pareto == True]
    df_biomass_pareto = df_yields.loc[is_pareto].sort_values('c_biomass_yield')
    df_biomass_pareto # efms in orange in figure below
    # compute convex hull
    from scipy.spatial import ConvexHull, convex_hull_plot_2d
    df_hull_all = ConvexHull(df_yields_filter[['c_biomass_yield', 'o2_biomass_yield']])
    df_hull_all.simplices
    sp = set()
    for simplex in df_hull_all.simplices:
        if (df_yields_filter.iloc[simplex[0]].name in df_biomass_pareto.index) and \
        (df_yields_filter.iloc[simplex[1]].name in df_biomass_pareto.index):
            sp.add(simplex[0])
            sp.add(simplex[1])
    sp = list(sp)
    df_hull = df_yields_filter.iloc[sp]
    df_hull = df_hull.sort_values('c_biomass_yield')
    df_hull
    df_pareto = df_hull.copy() # these are the efms we are interested in
    df_pareto # efms in red in figure above

    return df_pareto

