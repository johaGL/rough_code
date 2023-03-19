    parser.add_argument('--user_defined_path_to_prepared_tables', default=None,
                        help="by default,  'results/prepared_tables/' inside \
                             the folder defined as out_path in the configuration file")


#---------------


    lld = fg.read_clean_tables_names(f'{out_path}results/prepared_tables/TABLESNAMES.csv')

    
    # separately using fracContrib and abundance for the pca plots
    print(lld.keys())
    
    
    
    # seeing all tables in prepared : 
    
    for t in lld.keys():
        for co in compartments:
            papa = f'{out_path}results/prepared_tables/{lld[t]}--{co}--{suffix}.tsv'
            print(papa)
            print(os.path.exists(papa))
            
            
        print(fn)
        print(os.path.exists(fn))
        
        
        
        
 # test iris
    doiristest = True  # for debug
    if doiristest == True:
       
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    
