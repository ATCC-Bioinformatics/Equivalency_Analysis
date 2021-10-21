import pandas as pd
import seaborn as sns
import math
import matplotlib.pyplot as plt

df = pd.read_csv('ATCC_CheckM_stats_by_NCBI_ID.txt',sep='\t')

#contig plot
contdf = df[['Portal Contig Count', 'Genbank Contig Count']]
contdf['logGenbank'] = contdf['Genbank Contig Count'].apply(lambda x: math.log(x)) #logarithm used to prevent the genbank contig count max from overwhelming the color distribution
contdf['logPortal'] = contdf['Portal Contig Count'].apply(lambda x: math.log(x))
contiglog = contdf[['logPortal','logGenbank']]
sns.heatmap(contiglog,cmap='coolwarm').figure.savefig('logcont.png')

#completion plot
plt.clf() #remove older plot from memory so they don't overlap
plt.figure(figsize = (2.5,5))
compdf = df[['Portal Completion', 'GenBank Completion']]
sns.heatmap(compdf,vmin=85,vmax=100,).figure.savefig('completion.png') #a higher vmin makes this plot clearer. The lowest completion was 80.35.
#contamination plot
plt.clf()
plt.figure(figsize = (2.5,5))
contam = df[['Portal Contamination', 'GenBank Contamination']]
sns.heatmap(contam,cmap='coolwarm').figure.savefig('contamination.png')
#Remove entries where both have 0 contamination - attempt to make heatmap clearer. #Introduce robust coloring and remove y axis labels
contam_clean = contam[(contam['Portal Contamination'] != 0.000) & (contam['GenBank Contamination'] != 0.000)]
plt.clf()
plt.figure(figsize = (2.5,5))
sns.heatmap(contam_clean,cmap='coolwarm',robust=True,yticklabels=False).figure.savefig('contamination.nozeroes.robust.png')
#n50 plot
plt.clf()
plt.figure(figsize = (2.5,5))
n50 = df[['Portal N50', 'Genbank N50']]
sns.heatmap(n50,cmap='coolwarm').figure.savefig('n50.png')