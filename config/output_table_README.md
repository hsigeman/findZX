Output table README file
==========================

The output tables in this directory are used to generate output plots starting with numbers 1-4 (../plots/). The PDF version of each plot contain data paths to the tables that was used to create each one. 

All tables contain sex-specific genome coverage/heterozygosity values, as well as sex differences and ratios of these values. 



Column descriptions, tables with window values
----------------------------------------------

All tables with a window size (e.g. diffGenomeCoverage.mismatch.0.0.100000bp.out) contain the following columns:
1. chr (Chromosome/scaffold ID)
2. range (window number starting with 0)
3. heterogametic (mean value for heterogametic samples across all 5kb windows within the range)
4. homogametic (mean value for homogametic samples across all 5kb windows within the range)
5. ratio (mean ratio between the heterogametic and homogametic samples for all 5kb windows within the range)
6. ratio_scaled (same as column 4, but scaled by the mean)
7. diff (mean difference between the heterogametic and homogametic samples for all 5kb windows within the range)
8. diff_scaled (same as column 6, but scaled by the mean)
9. start (minimum chromosome/scaffold position within range)
10. end (maximum chromosome/scaffold position within range)



Column descriptions, tables with chromosome/scaffold values
-----------------------------------------------------------

All tables with per chromosome/scaffold values (e.g. diffGenomeCoverage.mismatch.0.0.chr.out) contain the following columns: 
1. chr (Chromosome/scaffold ID)
2. heterogametic (mean value for heterogametic samples across all 5kb windows within the range)
3. homogametic (mean value for homogametic samples across all 5kb windows within the range)
4. ratio (mean ratio between the heterogametic and homogametic samples for all 5kb windows within the range)
5. ratio_scaled (same as column 4, but scaled by the mean)
6. diff (mean difference between the heterogametic and homogametic samples for all 5kb windows within the range)
7. diff_scaled (same as column 6, but scaled by the mean)
8. length (approximate length of chromosome/scaffold; cumulative length of all windows belonging to this chromosome/scaffold)
