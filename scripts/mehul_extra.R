PlotTheme1 = theme_bw() +
            theme(axis.line = element_line(colour = "black"),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              strip.background = element_rect(colour="black", fill="white"),
              text = element_text(size = 12))

ntlist = c('-','A','G','C','T')

sample_levels = c('stock','passaged_stock','turbinate','sample','isolate') 

selectit = c('sample_name','type', 'gene_id','ntpos','varnt','aapos','varaa',
             'varcodon','varfreq', 'totalcount',
             'refnt','refaa', 'passaged_refnt','passaged_refaa',
             'mouse_id','sample_id','type2','sample_id2', 'mouse_id2'
            )

pull_cols = c('name','segment','ntpos','major','majorfreq',
                  'minor','minorfreq','binocheck','totalcount',
                  'aapos','majoraa','majorcodon',
                  'minoraa','minorcodon','gene_id') 


gene_list = c('5\'UTR',
                  'ORF1a',
                  'ORF1b',
                  'nsp1',
                  'nsp2',
                  'nsp3',
                  'nsp4',
                  'nsp5',
                  'nsp6',
                  'nsp7',
                  'nsp8',
                  'nsp9',
                  'nsp10',
                  'nsp11',
                  'nsp12a',
                  'nsp12b',
                  'nsp13',
                  'nsp14',
                  'nsp15',
                  'nsp16',
                  'S',
                  'ORF3a',
                  'E',
                  'M',
                  'ORF6',
                  'ORF7a',
                  'ORF7b',
                  'ORF8',
                  'N',
                  'ORF10',
                  '3\'UTR')


coding_list = c('ORF1a','ORF1b','S','ORF3a','E',
                'M','ORF6','ORF7a','ORF7b','ORF8','N','ORF10')
                
gene_ord_list = c('5\'UTR',
                  'ORF1a',
                  'ORF1b',
                  'nsp1',
                  'nsp2',
                  'nsp3',
                  'nsp4',
                  'nsp5',
                  'nsp6',
                  'nsp7',
                  'nsp8',
                  'nsp9',
                  'nsp10',
                  'nsp11',
                  'nsp12a',
                  'nsp12b',
                  'nsp13',
                  'nsp14',
                  'nsp15',
                  'nsp16',
                  'S',
                  'ORF3a',
                  'E',
                  'M',
                  'ORF6',
                  'ORF7a',
                  'ORF7b',
                  'ORF8',
                  'N',
                  'ORF10',
                  '3\'UTR',
                  'INTERGENIC')
nsps = c('nsp1',
                  'nsp2',
                  'nsp3',
                  'nsp4',
                  'nsp5',
                  'nsp6',
                  'nsp7',
                  'nsp8',
                  'nsp9',
                  'nsp10',
                  'nsp11',
                  'nsp12a',
                  'nsp12b',
                  'nsp13',
                  'nsp14',
                  'nsp15',
                  'nsp16'

)

orf_group = c('ORF1a','ORF1b')
spike = c("S")
end_groups = c('ORF3a',
                  'E',
                  'M',
                  'ORF6',
                  'ORF7a',
                  'ORF7b',
                  'ORF8',
                  'N',
                  'ORF10',
                  '3\'UTR')

aminoacids = c('G','A','L','M','F','W','K',
                'Q','E','S','P','V','I','C',
                'Y','H','R','N','D','T','*')

ntcolors = c('black','#fc8d62','#8d9fcb','#e88ac3','#a6d854')
names(ntcolors) = ntlist
nt_colScale_fill <- scale_fill_manual(name = "nt",values = ntcolors)
nt_colScale <- scale_colour_manual(name = "nt",values = ntcolors)

aacolors = as.vector(alphabet(21))
names(aacolors) = aminoacids
aa_colScale_fill <- scale_fill_manual(name = "aa",values = aacolors)
aa_colScale <- scale_colour_manual(name = "aa",values = aacolors)


gene_ord_list

col_list = c('black','#bdcee0','#dac8b9','#84bf3b',
                    '#66bace',
                    '#4372d6',
                    '#949cf3',
                    '#6233e3',
                    '#ba5fe5',
                    '#b03766',
                    '#db7032',
                    '#eece50',
                    '#cacaca',
            'black')
short_gene_list = c('5\'UTR','ORF1a','ORF1b','S','ORF3a','E','M','ORF6','ORF7a','ORF7b','ORF8','N','ORF10','3\'UTR')

names(col_list) = short_gene_list

gene_colScale_fill <- scale_fill_manual(name = "grp",values = col_list)

gene_colScale <- scale_colour_manual(name = "grp",values = col_list)


PullDels = function(files, del_list, del_minor, metadf, mergecol){
    del_positions = data.frame()
    MIN_del_positions = data.frame()
    
    
    for (f in files){
        df = read.csv(f, header = TRUE)
        df$name = gsub("-","_", df$name)

        temp = df %>%
                filter(ntpos %in% del_list$ntpos) %>%
            unique() # filter
        del_positions = rbind(del_positions, temp) # append

        temp = df %>%
                filter(ntpos %in% del_minor$ntpos) %>%
            unique()

        MIN_del_positions = rbind(MIN_del_positions, temp) # append
    }
    
    del_positions = del_positions %>% select(name, segment, ntpos,
                             major, majorfreq, minor, minorfreq, totalcount,
                             refnt,gene_id) %>%
    unique() 
    
    MIN_del_positions = MIN_del_positions %>% select(name, segment, ntpos,
                             major, majorfreq, minor, minorfreq, totalcount,
                             refnt,gene_id) %>%
                        unique() 
    
    del_positions = merge(del_positions, metadf, by.x = 'name', by.y = mergecol)
    MIN_del_positions = merge(MIN_del_positions, metadf, by.x = 'name', by.y = mergecol)
    
    return_list = list("del_positions" = del_positions, 
                  "MIN_del_positions" = MIN_del_positions)
    
    
    return(return_list)
    
}

plot_del = function(dels, orf_group, ntlist, major_cov){
    
    options(repr.plot.width = 17, repr.plot.height = 7)
    plotit = ggplot() + 

    
        geom_point(data = dels %>% 
                           filter(!gene_id %in% orf_group &
                                  major %in% ntlist &
                                  major != '-' &
                                  totalcount >= major_cov) %>% unique(), 
                   aes(x=factor(ntpos), y = sample_name, color = major),
                   shape = 'square', size =3) + 

        geom_point(data = dels %>% 
                           filter(!gene_id %in% orf_group &
                                  major %in% ntlist &
                                 passaged_stocknt != '-') %>% 
                           select(ntpos, gene_id, passaged_stocknt) %>% 
                           unique(), 
                   aes(x=factor(ntpos), y = "passaged_stocknt", color = passaged_stocknt),
                   shape = 'square', size =3) + 

        geom_point(data = dels %>% 
                           filter(!gene_id %in% 
                                  orf_group &
                                  passaged_stocknt == '-') %>% 
                           select(ntpos, gene_id, passaged_stocknt) %>% 
                           unique(), 
                   aes(x=factor(ntpos), y = "passaged_stocknt"),
                   shape = 'square', size =4, , color = 'black') + 

        geom_point(data = dels %>% 
                           filter(!gene_id %in% orf_group &
                                  major %in% ntlist &
                                 stocknt != '-') %>% 
                           select(ntpos, gene_id, stocknt) %>% 
                           unique(), 
                   aes(x=factor(ntpos), y = "stocknt", color = stocknt),
                   shape = 'square', size =3) + 

        geom_point(data = dels %>% 
                           filter(!gene_id %in% 
                                  orf_group &
                                  stocknt == '-') %>% 
                           select(ntpos, gene_id, stocknt) %>% 
                           unique(), 
                   aes(x=factor(ntpos), y = "stocknt"),
                   shape = 'square', size =4, , color = 'black') + 


        geom_point(data = dels %>% 
                           filter(!gene_id %in% orf_group &
                               major == "-" &
                                  totalcount >= major_cov) %>% 
                       unique(), 
                   aes(x=factor(ntpos), y = sample_name),
                   shape = 'square', size =4, color = 'black') + 

        PlotTheme1 +
        nt_colScale + 
        facet_grid(.~gene_id, scales = 'free', space = 'free')

    print(plotit)
    
}

plot_gene_del = function(dels, gene, samps_w_dels, orf_group, ntlist, major_cov){    
        p2 = ggplot() + 
            geom_point(data = dels %>% 
                               filter(gene_id == gene &
                                      sample_name %in% samps_w_dels &
                                      !gene_id %in% orf_group &
                                      major %in% ntlist &
                                      totalcount >= major_cov) %>% unique(), 
                       aes(x=factor(ntpos), y = sample_id2, color = major),
                       shape = 'square', size =3) + 


            geom_point(data = dels %>% 
                               filter(gene_id == gene &
                                      sample_name %in% samps_w_dels &
                                       major == '-' &
                                      totalcount >= major_cov &
                                      major != stocknt) %>% 
                           unique(), 
                       aes(x=factor(ntpos), y = sample_id2),
                       shape = 'square', size =4, color = 'black') + 

            PlotTheme1 +
            nt_colScale +
            scale_y_discrete(limits=rev)
    
    return(p2)
}

calcp = function(covdf, dp, genome_size = 29903){
    # covdf needs to be merged with metadata (have name and sample name)  and totalcount cols
    
    temp = covdf  %>%
        filter(totalcount > dp) %>%
        unique() %>%
        group_by(name, sample_name, rep) %>%
        #group_by(name.x, name.y, sample_name) %>%
        tally() %>%
        mutate(dp = dp,
              percent_cov = n/genome_size) 
    
    return(temp)
}

prep_stock = function(stock_files, repnum, major_cov = 10, minor_cov = 200, minor_freq=0.01,ntlist){
    stock_con = as.data.frame(seq(1, 29903, by = 1))
    colnames(stock_con) = c('ntpos')
    stock_cov = data.frame()
    genes_df = data.frame()
    stock_var = data.frame()  # only interested in minor vars potentially
    dp_cols = c()
    
    
    for (f in stock_files){

        df = read.csv(f, header = TRUE) 
        df = df %>%
              rowwise() %>% 
              mutate(name = gsub("-","_",name)) %>%
              ungroup()

        sample_name  = c(levels(factor(df$name)))[1]  # pull the sample name    

        # pull stock coverage information
        temp_cov = df %>% 
                select(name, ntpos, totalcount) %>%
                unique() %>%
                arrange(ntpos) # pull cov info across genome (remove gene id info)
        stock_cov = rbind(stock_cov, temp_cov)  ## append

        # pull consensus information
        temp_con = df %>% 
                select(ntpos, major, gene_id, aapos, majoraa, totalcount) %>%
                rowwise() %>%
                mutate(dp = ifelse(totalcount >= major_cov, 1, 0)) %>% # if it passes = 1, if not = 0
                select(-totalcount) %>%
                unique() %>%
                arrange(ntpos)   
        colnames(temp_con) = c("ntpos", glue("{sample_name}_{repnum}"), "gene_id", "aapos", glue("{sample_name}aa_{repnum}"), glue("{sample_name}dp_{repnum}"))

        # pull colnames for dp to compare to later
        dp_cols = c(dp_cols,  glue("{sample_name}dp_{repnum}"))

        if (f == stock_files[1]){
            # if the sample is the first in the list, append by ntpos
            stock_con = merge(stock_con, temp_con, by = c('ntpos'), all = TRUE)
        }else{
            # otherwise, merge by ntpos, gene id and aapos
            stock_con = merge(stock_con, temp_con, by = c('ntpos','gene_id','aapos'), all = TRUE)
        }


        # pull out gene info (for cov figures)
        gene_df = df %>%
                    select(name, gene_id, ntpos) %>%
                    unique()
        genes_df = rbind(genes_df, gene_df)

        # pull out minor variant information (not pulling for all because we want to compare)
        temp_var = df %>%
                    filter(minor %in% ntlist &
                           minorfreq >= minor_freq &
                           totalcount >= minor_cov
                          ) %>%
                    unique()

        stock_var = rbind(stock_var, temp_var)
    }

    stock_con = stock_con %>% arrange(ntpos, gene_id, aapos)
    stock_var = stock_var %>% arrange(ntpos, gene_id, aapos)

    # grab min and max of genes (for plotting)
    genes_df = genes_df %>%
            group_by(name, gene_id) %>%
            mutate(gene_start = min(ntpos),
                  gene_end = max(ntpos)) %>%
    ungroup() %>%
    unique()

       
    return_list = list("stock_con" = stock_con, "stock_cov" = stock_cov, "genes_df" = genes_df, "stock_var" = stock_var, "dp_cols" = dp_cols)
    return(return_list)
}


prep_positions = function(files, stock_con, major_cov, minor_cov, minor_freq,
                       ntlist, aminoacids, nsps){

    conlist = c()
    minlist = c()
    del_list = data.frame()
    del_minor = data.frame()

    for (f in files){
        df = read.csv(f, header = TRUE)
        df$name = gsub("-","_", df$name)

        # grab deletion positions across all samples
        # CONSENSUS deletions
        temp_del = df %>% filter(major == '-' &
                               totalcount >= major_cov) %>%
                select(ntpos, major) %>%
                unique()
        del_list = rbind(del_list, temp_del) %>% unique()

        # MINOR deletions:
        temp_mindel = df %>% filter(minor == '-' &
                                   totalcount >= minor_cov &
                                   minorfreq >= minor_freq) %>%
                select(ntpos, minor) %>%
                unique()

        del_minor = rbind(del_minor, temp_mindel) %>% unique()


        # add in the stock information to compare to: 
        df = merge(df, stock_con, 
                   by = c('ntpos', 'gene_id', 'aapos'), 
                   all = TRUE)


        # pull consensus info - compared to stock seq
        temp_con = df %>%
            filter(major != stocknt &
                   major %in% ntlist &
                   stocknt %in% ntlist &
                   passaged_stocknt %in% ntlist &
                   totalcount >= major_cov &
                  !gene_id %in% nsps &
            totalcount >= major_cov ) %>%
            unique()

        conlist = c(conlist, temp_con$ntpos)

        # pull consensus info - compared to stock seq
        temp_min = df %>%
            filter(minor %in% ntlist &
                   stocknt %in% ntlist &
                   passaged_stocknt %in% ntlist & 
                   totalcount >= minor_cov &
                   minorfreq >= minor_freq) %>%
            unique()

        minlist = c(minlist, temp_min$ntpos)


    }
    
    return_list = list("conlist"= conlist, "minlist" = minlist, "del_list" = del_list, "del_minor" = del_minor)
    return(return_list)

}


side_cov = function(cov_df, genes_df, minor_cov, major_cov, gene_ord_list, nsps, replicate_info){
    p1 = cov_df %>%
        ggplot(., aes(group = sample_name, x = ntpos, y = log10(totalcount))) +
            geom_line() +
            geom_hline(yintercept = log10(minor_cov), linetype = 2, color = 'black') +
            geom_hline(yintercept = log10(major_cov), linetype = 2, color = 'black') +
            PlotTheme1 +
            ggtitle(replicate_info) +
            geom_segment(data = genes_df %>%
                             select(-ntpos) %>% 
                             unique() %>%
                             filter(gene_id %in% gene_ord_list &
                                   !gene_id %in% nsps &
                                    sample_name %in% cov_df$sample_name) %>% unique(), 
                             aes(x= gene_start, xend = gene_end, 
                                             y=log10(100000), yend = log10(100000), 
                                             color = gene_id), linewidth = 3) + 
            facet_grid(.~sample_name, space = 'free', scales = 'free') +
            scale_color_manual(values = as.vector(polychrome(32))) 
    
    return(p1)
}

stacked_cov = function(cov_df, genes_df, minor_cov, major_cov, gene_ord_list, nsps, replicate_info){
    p2 = cov_df %>%
        filter() %>%
        ggplot(., aes(group = sample_name, x = ntpos, y = log10(totalcount))) +
            geom_line() +
            geom_hline(yintercept = log10(minor_cov), linetype = 2, color = 'black') +
            geom_hline(yintercept = log10(major_cov), linetype = 2, color = 'black') +
            PlotTheme1 +
            ggtitle(replicate_info) +
            geom_segment(data = genes_df %>%
                             select(-ntpos) %>%
                             unique() %>%
                             filter(gene_id %in% gene_ord_list &
                                   !gene_id %in% nsps &
                                    sample_name %in% cov_df$sample_name) %>% unique(), 
                             aes(x= gene_start, xend = gene_end, 
                                             y=log10(100000), yend = log10(100000), 
                                             color = gene_id), linewidth = 3) + 
            facet_grid(sample_name~., space = 'free', scales = 'free') + 
            scale_color_manual(values = as.vector(polychrome(32))) 

    return(p2)
}


prepStockMinors = function(df1, df2, minor_freq, minor_cov, ntlist){
 # a list of columns to pull from each rep var

    merge_var = c('sample_name','type','segment','ntpos','aapos','gene_id') # COLUMNS TO MERGE BY
    # not merging by actual major or minor - but will filter and require that they are the same below: 
    

    stock_var = merge(df1, df2, by = c(all_of(merge_var))) %>%
        filter(major.x == major.y & # major and minor nt in both reps have to be the same at each nt
               minor.x == minor.y &
               major.x %in% ntlist &
               major.y %in% ntlist &
               minor.x %in% ntlist &
               minor.y %in% ntlist
              ) %>% # filter out everything that isn't the same or doesn't pass cutoffs
        unique() %>%
        rowwise() %>%
        # since everything is the same average and pull out info
        mutate(major = major.x,
              minor = minor.x,
              majorfreq = (majorfreq.x + majorfreq.y)/2,
              minorfreq = (minorfreq.x + minorfreq.y)/2,
              minorcodon = minorcodon.x,
              minoraa = minoraa.x,
              totalcount = (totalcount.x + totalcount.y)/2, 
              binocheck = ifelse(binocheck.x == 'True' | binocheck.y == 'True', 'PASS', 'FAIL')) %>%
        unique()
    
    return(stock_var)
}

summarizeMinors = function(df, orf_group, minor_freq, minor_cov){
    return_df = df %>%    
                    filter(!gene_id %in% orf_group &
                          minorfreq >= minor_freq &
                          totalcount >= minor_cov &
                          binocheck == 'PASS'
                          ) %>% 
                unique()
    
    return_df = return_df %>% 
                select(ntpos, major, minor, minorfreq, totalcount, gene_id) %>% 
                group_by(ntpos, minor, major, gene_id) %>%
                mutate(
                    mean_minorfreq = mean(minorfreq),
                    sd_minorfreq = sd(minorfreq),
                    min_minorfreq = min(minorfreq),
                    max_minorfreq = max(minorfreq),
                    avg_dp = mean(totalcount),
                    sd_dp = sd(totalcount),
                    sample_number = n()
                    ) %>%
                select(-minorfreq, -totalcount) %>%
                unique() %>%
                arrange(sample_number) %>%
                data.frame()
    
    return(return_df)

    
}


prep_var_dfs = function(files, stock_con, major_cov, minor_cov, minor_freq, conlist, stock_minor_summary, pass_stock_minor_summary, minlist, ntlist){
    consensus_df = data.frame()
    conlist = levels(factor(conlist))

    minor_df = data.frame()
    minlist = levels(factor(minlist))

    for (f in files){
        df = read.csv(f, header = TRUE)
        df$name = gsub("-","_", df$name)

        df = merge(df, stock_con, 
                   by = c('ntpos','gene_id','aapos'), 
                   all = TRUE)

        # pull consensus info - compared to stock seq
        temp_con = df %>%
            filter(totalcount >= major_cov &
                   ntpos %in% conlist |
                totalcount >= major_cov &
                   ntpos %in% stock_minor_summary$ntpos |
                totalcount >= major_cov &
                   ntpos %in% pass_stock_minor_summary$ntpos) %>% ## updated 06.23.2023 adding the stock minor positions to check
            unique()
        consensus_df = rbind(consensus_df, temp_con)


         # pull consensus info - compared to stock seq
        temp_min = df %>%
            filter(ntpos %in% minlist | 
                   ntpos %in% conlist | 
                   ntpos %in% stock_minor_summary$ntpos |
                   ntpos %in% pass_stock_minor_summary$ntpos) %>%
            unique()
        minor_df = rbind(minor_df, temp_con)

    }
    
    # CONSENSUS CHANGES ARE IN COMPARISON TO THE STOCK 
    consensus_df = consensus_df %>% 
        rowwise() %>%
        mutate(
               vartype = 'major',
               major = ifelse(major %in% ntlist &
                             totalcount >= major_cov, major, NA),
               majorfreq = ifelse(major %in% ntlist &
                                 totalcount >= major_cov, majorfreq, NA)) %>%
        
        select(name, major,ntpos, totalcount, majorfreq, 
               aapos, gene_id, majoraa, majorcodon, stocknt, stockaa, passaged_stocknt, passaged_stockaa, vartype) %>%
        unique() %>% 
        ungroup() %>%
        rename("varnt" = "major",
                "varfreq" = "majorfreq",
                "varaa" = "majoraa",
                "varcodon" = "majorcodon",
                "refnt" = "stocknt", 
                 "refaa" = "stockaa",
              "passaged_refnt" = "passaged_stocknt",
              "passaged_refaa" = "passaged_stockaa") 
    
    # MINORS ARE IN COMPARISON TO THE MAJOR INFORMATION IN THE SAMPLE
    minor_df = minor_df %>% 
        rowwise() %>%
        mutate(
              vartype = 'minor',
              minor = ifelse(minor %in% ntlist & 
                            minorfreq >= minor_freq &
                            totalcount >= minor_cov, 
                             minor, NA),
              minorfreq = ifelse(minor %in% ntlist & 
                            minorfreq >= minor_freq &
                            totalcount >= minor_cov, minorfreq, NA)) %>%
        select(name, minor, ntpos, totalcount, minorfreq, 
               aapos, gene_id, minoraa, minorcodon,stocknt, stockaa, passaged_stocknt, passaged_stockaa, vartype) %>%
        unique() %>% 
        ungroup() %>%
        rename("varnt" = "minor",
                "varfreq" = "minorfreq",
                "varaa" = "minoraa",
                "varcodon" = 'minorcodon',
               "refnt" = "stocknt", 
                 "refaa" = "stockaa",
              "passaged_refnt" = "passaged_stocknt",
              "passaged_refaa" = "passaged_stockaa") 
    
    return_list = list("consensus_df" = consensus_df, 
                  "minor_df" = minor_df)
    
    return(return_list)
}

# addedon 06.21.2023 to adjust for everything
cov_sample_prep = function(files, meta){
    cov_df = data.frame()
    genes_df = data.frame()


    for (f in files){
        df = read.csv(f, header = TRUE)  %>%
              rowwise() %>% 
              mutate(name = gsub("-","_",name)) %>%
              ungroup()
        
        temp_cov = df %>% 
                select(name, ntpos, totalcount) %>%
                unique() %>%
                arrange(ntpos)  # temp coverage info for ach sample
        cov_df = rbind(cov_df, temp_cov)  ## append


        # pull out gene info (for figures)
        gene_df = df %>%
                    select(name, gene_id, ntpos) %>%
                    unique()
        genes_df = rbind(genes_df, gene_df)

    }

    genes_df = genes_df %>%
                group_by(name, gene_id) %>%
                mutate(gene_start = min(ntpos),
                      gene_end = max(ntpos)) %>%
        ungroup() %>%
        unique()
    
    return_list = list("cov_df" = cov_df, "genes_df" = genes_df)

    return(return_list)
}

sample_labels = c('1_input','2_input','3_input','4_input',
                  '1_input_vero','2_input_vero','3_input_vero',
                  '1_input-C57/BL6','2_input-C57/BL6','3_input-C57/BL6','4_input-C57/BL6','5_input-C57/BL6',
                  '21_sample','22_sample','23_sample','24_sample','25_sample',
                  '26_sample','27_sample','28_sample','29_sample','30_sample',
                  '24_isolate','26_isolate','27_isolate','29_isolate','30_isolate')



name_labels = c('B1_351_INPUT_1','B1_351_INPUT_2','B1_351_INPUT_3','B1_351_INPUT_4',
                'B_1_351_B1','B_1_351_B2','B_1_351_B3',
                'MK_4A_1_RNA_turbinate','MK_4A_2_RNA_turbinate','MK_4A_3_RNA_turbinate','MK_4A_4_RNA_turbinate','MK_4A_5_RNA_turbinate',
                'KAF_109_21','KAF_109_22','KAF_109_23','KAF_109_24','KAF_109_25',
                'KAF_109_26','KAF_109_27','KAF_109_28','KAF_109_29','KAF_109_30',
                'isolate_24','isolate_26','isolate_27','isolate_29','isolate_30')


name_colors = c('#66a5c8','#4d94c0','#3285b7','#0067A5',
                  '#e0d0e5','#c497cf','#ac90c3',
                  '#ffa5d2','#ff96cb','#ff87c3','#ff78bc','#ff69b4',
                '#848482','#222222','#F3C300','#f94040','#f38400',
                '#00c000','#0000ff','#00998f','#ad07e3','#55a0fb',
                '#f94040','#00c000','#0000ff','#ad07e3','#55a0fb')

names(name_colors) = name_labels
name_colScale_fill <- scale_fill_manual(name = "id",values = name_colors)
name_colScale <- scale_colour_manual(name = "id",values = name_colors)

sample_colors = name_colors
names(sample_colors) = sample_labels
sample_colScale_fill <- scale_fill_manual(name = "id",values = sample_colors)
sample_colScale <- scale_colour_manual(name = "id",values = sample_colors)

id_labels = c('21','22','23','24','25',
              '26','27','28','29','30',
              'input-vero','input','input-c57/bl6')


id_colors = c('#848482','#222222','#F3C300','#f94040','#f38400',
              '#00c000','#0000ff','#00998f','#ad07e3','#55a0fb',
              '#c497cf','#0067A5','#ff69b4')


names(id_colors) = id_labels
id_colScale_fill <- scale_fill_manual(name = "id",values = id_colors)
id_colScale <- scale_colour_manual(name = "id",values = id_colors)


type_names = c('input','input_vero','input-C57/BL6','isolate','sample')
shape_numbers = c(18, 15, 17, 4, 16)
names(shape_numbers) = type_names
type_shapes = scale_shape_manual(name = 'type',values = shape_numbers)

PlotIsolates = function(isolate_number, DF){
    isolate_df = DF %>% filter(mouse_id == isolate_number & varnt != refnt & vartype == 'major') %>%
                    group_by(ntpos) %>%
                    add_tally(name='isolate_count') %>%
                    ungroup()
    
    plotit = ggplot() + 
        geom_vline(data = isolate_df %>% 
                   filter(isolate_count == 2), 
                   aes(xintercept = ntpos), color = 'black', linetype = 2 ) + 
        
        # in the original mouse as a consensus change but not in the isolate-vero sample
        geom_vline(data = isolate_df %>% 
                   filter(isolate_count == 1 & type == 'sample'), 
                   aes(xintercept = ntpos), color = 'gray', linetype = 2 ) + 
        
        # not in the original mouse as a consensus change but in the isolate-vero sample as a consensus change
        geom_vline(data = isolate_df %>% 
                   filter(isolate_count == 1 & type == 'isolate'), 
                   aes(xintercept = ntpos), color = 'blue', linetype = 2 ) + 

        geom_text(data = isolate_df %>% filter(isolate_count == 2 &
                                              varaa %in% aminoacids &
                                              refaa %in% aminoacids) %>%
                            select(ntpos, gene_id, var_info) %>% 
                            unique(),
                  aes(x = ntpos, y = 1.15, label = glue("{gene_id}: {var_info}")), color = 'black', angle = 90, size = 3) + 

        geom_text(data = isolate_df %>% filter(isolate_count == 1 & type == 'sample' &
                                              varaa %in% aminoacids &
                                              refaa %in% aminoacids) %>%
                            select(ntpos, gene_id, var_info) %>% 
                            unique(),
                  aes(x = ntpos, y = 1.15, label = glue("{gene_id}: {var_info}")),  color = 'gray', angle = 90, size = 3) + 

        geom_text(data = isolate_df %>% filter(isolate_count == 1 & type == 'isolate' &
                                              varaa %in% aminoacids &
                                              refaa %in% aminoacids) %>%
                            select(ntpos, gene_id, var_info) %>% 
                            unique(),
                  aes(x = ntpos, y = 1.15, 
                      label = glue("{gene_id}: {var_info}")), color = 'blue',angle = 90, size = 3) + 
    
    
    # nonsyn or del variants
    geom_text(data = isolate_df %>% filter(isolate_count == 2) %>% 
                                     filter(!varaa %in% aminoacids | 
                                              !refaa %in% aminoacids) %>%
                            select(ntpos, gene_id, var_info_nt) %>% 
                            unique(),
                  aes(x = ntpos, y = 1.15, label = glue("{gene_id}: {var_info_nt}")), color = 'black', angle = 90, size = 3) + 

        geom_text(data = isolate_df %>% filter(isolate_count == 1 & type == 'sample') %>% 
                                     filter(!varaa %in% aminoacids | 
                                              !refaa %in% aminoacids) %>%
                            select(ntpos, gene_id, var_info_nt) %>% 
                            unique(),
                  aes(x = ntpos, y = 1.15, label = glue("{gene_id}: {var_info_nt}")),  color = 'gray', angle = 90, size = 3) + 

        geom_text(data = isolate_df %>% filter(isolate_count == 1 & type == 'isolate') %>% 
                                     filter(!varaa %in% aminoacids | 
                                              !refaa %in% aminoacids) %>%
                            select(ntpos, gene_id, var_info_nt) %>% 
                            unique(),
                  aes(x = ntpos, y = 1.15, 
                      label = glue("{gene_id}: {var_info_nt}")), color = 'blue',angle = 90, size = 3) + 
    
    
    
    
        PlotTheme1 + 
        geom_gene_arrow(data = gene_info %>% filter(!NAME %in% nsps) , 
                   aes(xmin = START, xmax = END, y = 1.1, fill = NAME)) +
        gene_colScale_fill + 
        ggtitle(glue("Isolate: {isolate_number}, black = both, gray=mouse only, blue=isolate only")) + 
        PlotTheme1 +
        theme(legend.position="bottom") +
        ylim(1, 1.25)
    
    return(plotit)
}
