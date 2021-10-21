import pandas as pd

df_pub = pd.read_csv('Portal_public_spreadsheet_public.txt',sep='\t')
df_pub = df_pub[df_pub['ATCC Assembled  yes/no'].str.lower() == 'yes'] #get only those records which have an atcc assembly
df_pub['ATCC Catalog'] = df_pub['ATCC Catalog'].apply(lambda x: x.upper()) #set all catalogIDs to uppercase #count = 577
df_port = pd.read_csv('Portal_public_spreadsheet_portal.txt',sep='\t')
df_port['Base Catalog Number'] = df_port['Base Catalog Number'].apply(lambda x: x.upper()) #set all catalogIDs to uppercase

portal_catalog = list(set(list((df_port['Base Catalog Number'])))) #get list of pcatalog IDs on the portal. MAY BE REDUNDANT

#plan: return ratio values of length, illumina n50, contig count, read depth, gc content - gc content not currently in public dataset

#int(df_pub[df_pub['ATCC Catalog'] == catalog]['Total Length']) / int(df_port[df_port['Base Catalog Number'] == catalog]['Total Length']) #need to remove redundancy problem
#ideas:
#separate into assembly level data frames - still has redundancies by catalog id.
#return average ratio values when catalog ID appears more than once
#Both of the above so values returned are collected by assembly level. Would need to prevent lower levels from repeating higher level assemblies per catalogID
#third option will be addressed below

#set public assembly dataframes
df_pub_complete = df_pub[(df_pub['Assembly Level'] == 'Chromosome') | (df_pub['Assembly Level'] == 'Complete Genome')] #count=234
pub_complete = [x for x in portal_catalog if x in list(df_pub_complete['ATCC Catalog'])] #get "complete" catalog IDs from public assemblies
df_pub_scaffold = df_pub[~df_pub['ATCC Catalog'].isin(pub_complete)] #remove those catalogIDs that are considered complete
df_pub_scaffold = df_pub_scaffold[df_pub_scaffold['Assembly Level'] == 'Scaffold'] #get down to just scaffolds
pub_scaffold = [x for x in portal_catalog if x in list(df_pub_scaffold['ATCC Catalog'])] #get list of scaffold catalogIDs #count = 107
df_pub_contig = df_pub[~(df_pub['ATCC Catalog'].isin(pub_complete)) & ~(df_pub['ATCC Catalog'].isin(pub_scaffold))] #remove catalogIDs that are scaffold or complete #really doesn't need to go back to the top-level df, but it works
df_pub_contig = df_pub_contig[df_pub_contig['Assembly Level'] == 'Contig'] #just contigs #should be unnecessary at this stage #count = 81
#totals of dataframes don't add up to top level dataframe count because of assemblies that exist at multiple stages. This approach should only compare the highest quality public assembly to ATCC's.
df_pub_comp_mean = df_pub_complete.groupby('ATCC Catalog').mean().reset_index() #summarize values by average
df_pub_comp_min = df_pub_complete.groupby('ATCC Catalog').min().reset_index() #summarize values by minimum (appropriate for contig counts)
df_pub_comp_max = df_pub_complete.groupby('ATCC Catalog').max().reset_index() #summarize values by maximum (appropriate for N50, perhaps)

df_pub_scaffold_mean = df_pub_scaffold.groupby('ATCC Catalog').mean().reset_index() #summarize values by average
df_pub_scaffold_min = df_pub_scaffold.groupby('ATCC Catalog').min().reset_index() #summarize values by minimum (appropriate for contig counts)
df_pub_scaffold_max = df_pub_scaffold.groupby('ATCC Catalog').max().reset_index() #summarize values by maximum (appropriate for N50, perhaps)

df_pub_contig_mean = df_pub_contig.groupby('ATCC Catalog').mean().reset_index() #summarize values by average
df_pub_contig_min = df_pub_contig.groupby('ATCC Catalog').min().reset_index() #summarize values by minimum (appropriate for contig counts)
df_pub_contig_max = df_pub_contig.groupby('ATCC Catalog').max().reset_index() #summarize values by maximum (appropriate for N50, perhaps)


#set portal assembly dataframes
df_port_complete = df_port[df_port['Assembly Rank'] == 1] #count=234
if len(df_port_complete) == len(portal_catalog):
    print('Evaluating only top ranked assemblies')
    if len(list(set(list(df_port_complete['Base Catalog Number'])))) == len(portal_catalog):
        print('No Repeated Catalog IDs will be evaluated')
    else:
        print('There are Catalog IDs are present more than once in the top ranked assemblies, which will be averaged')
else:
    port_cat_2 = list(df_port[df_port['Assembly Rank'] == 2]['Base Catalog Number']) #get all second rank assemblies to check that each assembly has a top ranking
    if len([x for x in port_cat_2 if x not in list(df_port[df_port['Assembly Rank'] == 1]['Base Catalog Number'])]) == 0:
        print('No second rank assemblies without a first rank assembly in the portal assembly list')

#comparisons
#'complete' public assemblies
#genome length #compared by average length #comparison metric could be command line argument and consistent across all fields.
port_length = df_port_complete[['Base Catalog Number', 'Total Length', 'Product Collection']]
pub_length = df_pub_complete[['ATCC Catalog', 'Total Length', 'Assembly Accession', 'Organism']]
port_length.columns = ['Catalog ID', 'Total Length', 'Product Collection']
pub_length.columns = ['Catalog ID', 'Total Length', 'Assembly Accession', 'Organism'] #need same column names to merge dataframes
df_len = pd.merge(port_length, pub_length, on=['Catalog ID'], suffixes=('_portal','_public'))
df_len['Length Ratio'] = df_len['Total Length_portal'] / df_len['Total Length_public']
df_len['Length Difference'] = df_len['Total Length_portal'] - df_len['Total Length_public']
df_len['Public Assembly Level'] = 'Complete'
df_len.to_csv('Assembly_Complete_length_no_aggregate.comparisons.txt',sep='\t',index=False)

#N50 #compared by max N50
port_N50= df_port_complete[['Base Catalog Number', 'Filtered N50', 'Product Collection']]
pub_N50 = df_pub_complete[['ATCC Catalog', 'Scaffold N50', 'Assembly Accession', 'Organism']]
port_N50.columns = ['Catalog ID', 'N50', 'Product Collection']
pub_N50.columns = ['Catalog ID', 'N50', 'Assembly Accession', 'Organism']

df_N50 = pd.merge(port_N50, pub_N50, on=['Catalog ID'], suffixes=('_portal','_public'))
df_N50['N50 Ratio'] = df_N50['N50_portal'] / df_N50['N50_public']
df_N50['N50 Difference'] = df_N50['N50_portal'] - df_N50['N50_public']
df_N50['Public Assembly Level'] = 'Complete'
df_N50.to_csv('Assembly_Complete_N50_no_aggregate.comparisons.txt',sep='\t',index=False)

#contig/replicon count #compared by min
port_contig= df_port_complete[['Base Catalog Number', 'Total Contigs', 'Product Collection']]
pub_contig = df_pub_complete[['ATCC Catalog', 'Scaffold Count', 'Assembly Accession', 'Organism']]
port_contig.columns = ['Catalog ID', 'Contig Count', 'Product Collection']
pub_contig.columns = ['Catalog ID', 'Contig Count', 'Assembly Accession', 'Organism']

df_contig = pd.merge(port_contig, pub_contig, on=['Catalog ID'], suffixes=('_portal','_public'))
df_contig['Contig Count Ratio'] = df_contig['Contig Count_portal'] / df_contig['Contig Count_public']
df_contig['Contig Count Difference'] = df_contig['Contig Count_portal'] - df_contig['Contig Count_public']
df_contig['Public Assembly Level'] = 'Complete'
df_contig.to_csv('Assembly_Complete_Contig_no_aggregate.comparisons.txt',sep='\t',index=False)

#GC Content #compared by mean
port_gc= df_port_complete[['Base Catalog Number', 'GC Content', 'Product Collection']]
pub_gc = df_pub_complete[['ATCC Catalog', 'GC Content', 'Assembly Accession', 'Organism']]
port_gc.columns = ['Catalog ID', 'GC', 'Product Collection']
port_gc['GC'] = port_gc['GC'] * 100
pub_gc.columns = ['Catalog ID', 'GC', 'Assembly Accession', 'Organism']

df_gc = pd.merge(port_gc, pub_gc, on=['Catalog ID'], suffixes=('_portal','_public'))
df_gc['GC Ratio'] = df_gc['GC_portal'] / df_gc['GC_public']
df_gc['GC Difference'] = df_gc['GC_portal'] - df_gc['GC_public']
df_gc['Public Assembly Level'] = 'Complete'
df_gc.to_csv('Assembly_Complete_GC_no_aggregate.comparisons.txt',sep='\t',index=False)

#'scaffold' public assemblies #note: variables are reused because we're writing to files as we go
#genome length #compared by average length #comparison metric could be command line argument and consistent across all fields.
#port_length = df_port_complete[['Base Catalog Number', 'Total Length']] #already done, doesn't need to be repeated
pub_length = df_pub_scaffold[['ATCC Catalog', 'Total Length', 'Assembly Accession', 'Organism']]
#port_length.columns = ['Catalog ID', 'Total Length'] #already done, doesn't need to be repeated
pub_length.columns = ['Catalog ID', 'Total Length', 'Assembly Accession', 'Organism'] #need same column names to merge dataframes
df_len = pd.merge(port_length, pub_length, on=['Catalog ID'], suffixes=('_portal','_public'))
df_len['Length Ratio'] = df_len['Total Length_portal'] / df_len['Total Length_public']
df_len['Length Difference'] = df_len['Total Length_portal'] - df_len['Total Length_public']
df_len['Public Assembly Level'] = 'Scaffold'
df_len.to_csv('Assembly_Scaffold_length_no_aggregate.comparisons.txt',sep='\t',index=False)

#N50 #compared by max N50
#port_N50= df_port_complete[['Base Catalog Number', 'Filtered N50']]
pub_N50 = df_pub_scaffold[['ATCC Catalog', 'Scaffold N50', 'Assembly Accession', 'Organism']]
#port_N50.columns = ['Catalog ID', 'N50']
pub_N50.columns = ['Catalog ID', 'N50', 'Assembly Accession', 'Organism']

df_N50 = pd.merge(port_N50, pub_N50, on=['Catalog ID'], suffixes=('_portal','_public'))
df_N50['N50 Ratio'] = df_N50['N50_portal'] / df_N50['N50_public']
df_N50['N50 Difference'] = df_N50['N50_portal'] - df_N50['N50_public']
df_N50['Public Assembly Level'] = 'Scaffold'
df_N50.to_csv('Assembly_Scaffold_N50_no_aggregate.comparisons.txt',sep='\t',index=False)

#contig/replicon count #compared by min
#port_contig= df_port_complete[['Base Catalog Number', 'Total Contigs']]
pub_contig = df_pub_scaffold[['ATCC Catalog', 'Scaffold Count', 'Assembly Accession', 'Organism']]
#port_contig.columns = ['Catalog ID', 'Contig Count']
pub_contig.columns = ['Catalog ID', 'Contig Count', 'Assembly Accession', 'Organism']

df_contig = pd.merge(port_contig, pub_contig, on=['Catalog ID'], suffixes=('_portal','_public'))
df_contig['Contig Count Ratio'] = df_contig['Contig Count_portal'] / df_contig['Contig Count_public']
df_contig['Contig Count Difference'] = df_contig['Contig Count_portal'] - df_contig['Contig Count_public']
df_contig['Public Assembly Level'] = 'Scaffold'
df_contig.to_csv('Assembly_Scaffold_Contig_no_aggregate.comparisons.txt',sep='\t',index=False)

#GC Content #compared by mean
#port_gc= df_port_complete[['Base Catalog Number', 'GC Content']]
pub_gc = df_pub_scaffold[['ATCC Catalog', 'GC Content','Assembly Accession', 'Organism']]
#port_gc.columns = ['Catalog ID', 'GC']
#port_gc['GC'] = port_gc['GC'] * 100
pub_gc.columns = ['Catalog ID', 'GC', 'Assembly Accession', 'Organism']

df_gc = pd.merge(port_gc, pub_gc, on=['Catalog ID'], suffixes=('_portal','_public'))
df_gc['GC Ratio'] = df_gc['GC_portal'] / df_gc['GC_public']
df_gc['GC Difference'] = df_gc['GC_portal'] - df_gc['GC_public']
df_gc['Public Assembly Level'] = 'Scaffold'
df_gc.to_csv('Assembly_Scaffold_GC_no_aggregate.comparisons.txt',sep='\t',index=False)

#'contig' public assemblies
#genome length #compared by average length #comparison metric could be command line argument and consistent across all fields.
#port_length = df_port_complete[['Base Catalog Number', 'Total Length']]
pub_length = df_pub_contig[['ATCC Catalog', 'Total Length', 'Assembly Accession', 'Organism']]
#port_length.columns = ['Catalog ID', 'Total Length']
pub_length.columns = ['Catalog ID', 'Total Length', 'Assembly Accession', 'Organism'] #need same column names to merge dataframes
df_len = pd.merge(port_length, pub_length, on=['Catalog ID'], suffixes=('_portal','_public'))
df_len['Length Ratio'] = df_len['Total Length_portal'] / df_len['Total Length_public']
df_len['Length Difference'] = df_len['Total Length_portal'] - df_len['Total Length_public']
df_len['Public Assembly Level'] = 'Contig'
df_len.to_csv('Assembly_Contig_length_no_aggregate.comparisons.txt',sep='\t',index=False)

#N50 #compared by max N50
#port_N50= df_port_complete[['Base Catalog Number', 'Filtered N50']]
pub_N50 = df_pub_contig[['ATCC Catalog', 'Scaffold N50', 'Assembly Accession', 'Organism']]
#port_N50.columns = ['Catalog ID', 'N50']
pub_N50.columns = ['Catalog ID', 'N50', 'Assembly Accession', 'Organism']

df_N50 = pd.merge(port_N50, pub_N50, on=['Catalog ID'], suffixes=('_portal','_public'))
df_N50['N50 Ratio'] = df_N50['N50_portal'] / df_N50['N50_public']
df_N50['N50 Difference'] = df_N50['N50_portal'] - df_N50['N50_public']
df_N50['Public Assembly Level'] = 'Contig'
df_N50.to_csv('Assembly_Contig_N50_no_aggregate.comparisons.txt',sep='\t',index=False)

#contig/replicon count #compared by min
#port_contig= df_port_complete[['Base Catalog Number', 'Total Contigs']]
pub_contig = df_pub_contig[['ATCC Catalog', 'Scaffold Count', 'Assembly Accession','Organism']]
#port_contig.columns = ['Catalog ID', 'Contig Count']
pub_contig.columns = ['Catalog ID', 'Contig Count', 'Assembly Accession', 'Organism']

df_contig = pd.merge(port_contig, pub_contig, on=['Catalog ID'], suffixes=('_portal','_public'))
df_contig['Contig Count Ratio'] = df_contig['Contig Count_portal'] / df_contig['Contig Count_public']
df_contig['Contig Count Difference'] = df_contig['Contig Count_portal'] - df_contig['Contig Count_public']
df_contig['Public Assembly Level'] = 'Contig'
df_contig.to_csv('Assembly_Contig_Contig_no_aggregate.comparisons.txt',sep='\t',index=False)

#GC Content #compared by mean
#port_gc= df_port_complete[['Base Catalog Number', 'GC Content']]
pub_gc = df_pub_contig[['ATCC Catalog', 'GC Content', 'Assembly Accession', 'Organism']]
#port_gc.columns = ['Catalog ID', 'GC']
#port_gc['GC'] = port_gc['GC'] * 100
pub_gc.columns = ['Catalog ID', 'GC', 'Assembly Accession', 'Organism']

df_gc = pd.merge(port_gc, pub_gc, on=['Catalog ID'], suffixes=('_portal','_public'))
df_gc['GC Ratio'] = df_gc['GC_portal'] / df_gc['GC_public']
df_gc['GC Difference'] = df_gc['GC_portal'] - df_gc['GC_public']
df_gc['Public Assembly Level'] = 'Contig'
df_gc.to_csv('Assembly_Contig_GC_no_aggregate.comparisons.txt',sep='\t',index=False)