library(sjPlot)
library(ggplot2)
library(cowplot)

set.seed(999)

source("workflow/scripts/functions.R")
args <- commandArgs(trailingOnly = TRUE)

file1 = args[1]
file2 = args[2]
file3 = args[3]
filesnp = args[4]
plot_out = args[5]
table_out = args[6]
chr_file = args[7]
highlight_file = args[8]
ED1 = args[9]
ED2 = args[10]
ED3 = args[11]
WINDOW = args[12]
CHR_NR = args[13]

ED1 = gsub("\\.", "-", ED1)
ED2 = gsub("\\.", "-", ED2)
ED3 = gsub("\\.", "-", ED3)
ED1 = gsub("$", " mismatches", ED1)
ED2 = gsub("$", " mismatches", ED2)

ED1 = gsub("0-0", "0", ED1)
ED1 = gsub("0-", "<= ", ED1)
ED2 = gsub("0-", "<= ", ED2)

file1_base = gsub(".out$", "", file1)
file2_base = gsub(".out$", "", file2)
file3_base = gsub(".out$", "", file3)
filesnp_base = gsub(".out$", "", filesnp)

plot_out_base = gsub(".pdf$", "", plot_out)
table_out_base = gsub(".html$", "", table_out)

################################################################################
################################# READ FILES ###################################
################################################################################

ED1 = gsub("\\.", "-", ED1)
ED2 = gsub("\\.", "-", ED2)
ED3 = gsub("\\.", "-", ED3)

cov_1_table <- read.table( file1, 
                           header=TRUE, 
                           fill=TRUE, 
                           stringsAsFactor=FALSE)

cov_2_table <- read.table( file2, 
                           header=TRUE, 
                           fill=TRUE, 
                           stringsAsFactor=FALSE)

cov_3_table <- read.table( file3, 
                           header=TRUE, 
                           fill=TRUE, 
                           stringsAsFactor=FALSE)

snp_table   <- read.table( filesnp, 
                           header=TRUE, 
                           fill=TRUE, 
                           stringsAsFactor=FALSE)

if ( dim( cov_1_table )[1] == 0) {
  
  print("Warning: No chromosome/scaffold longer than the specified window size! Check output based 
        on chromosome instead.")
  
} else {
  
  if ( file.exists( chr_file ) ) {
    
    #chromosome <- read.csv(chr_file, 
    #                       header = FALSE, 
    #                       sep = ",")
    chromosome <- read.delim(chr_file, header = FALSE)
    chromosome$V1 <- trimws(chromosome$V1, which = c("both", "left", "right"), whitespace = "[ \t\r\n]")
    chromosome <- chromosome$V1
    cov_1_table <- cov_1_table[ cov_1_table$chr %in% chromosome, ]
    
  }
  
  nr_factors <- c( length( unique( cov_1_table$chr ) ),
                   length( unique( cov_2_table$chr ) ),
                   length( unique( cov_3_table$chr ) ),
                   length( unique( snp_table$chr )) )
  
  # Makes sure that all tables have the same chromosomes
  while ( min( nr_factors ) != max( nr_factors ) ) {
    
    nr_factors_min <- which.min(nr_factors)
    
    if ( nr_factors_min == 1 ) {
      
      cov_2_table <- cov_2_table[ cov_2_table$chr %in% unique(cov_1_table$chr), ]
      cov_3_table <- cov_3_table[ cov_3_table$chr %in% unique(cov_1_table$chr), ]
      snp_table <- snp_table[ snp_table$chr %in% unique(cov_1_table$chr), ]
      
    } else if ( nr_factors_min == 2 ) {
      
      cov_1_table <- cov_1_table[ cov_1_table$chr %in% unique(cov_2_table$chr), ]
      cov_3_table <- cov_3_table[ cov_3_table$chr %in% unique(cov_2_table$chr), ]
      snp_table <- snp_table[ snp_table$chr %in% unique(cov_2_table$chr), ]
      
    } else if ( nr_factors_min == 3 ) {
      
      cov_1_table <- cov_1_table[ cov_1_table$chr %in% unique(cov_3_table$chr), ]
      cov_2_table <- cov_2_table[ cov_2_table$chr %in% unique(cov_3_table$chr), ]
      snp_table <- snp_table[ snp_table$chr %in% unique(cov_3_table$chr), ]
      
    } else if ( nr_factors_min == 4 ) {
      
      cov_1_table <- cov_1_table[ cov_1_table$chr %in% unique(snp_table$chr), ]
      cov_2_table <- cov_2_table[ cov_2_table$chr %in% unique(snp_table$chr), ]
      cov_3_table <- cov_3_table[ cov_3_table$chr %in% unique(snp_table$chr), ]
      
    }
    
    nr_factors <- c( length(unique(cov_1_table$chr)),
                     length(unique(cov_2_table$chr)),
                     length(unique(cov_3_table$chr)),
                     length(unique(snp_table$chr)))
    
  }
  
  ################################################################################
  ############################# STATISTIC TEST ################################
  ################################################################################
  
  cov_1_table$diff_Zscaled <- (cov_1_table$diff - mean(cov_1_table$diff)) / sd(cov_1_table$diff)
  fit_cov1 <- lm(diff_Zscaled ~ chr -1, data = cov_1_table)
  
  cov_2_table$diff_Zscaled <- (cov_2_table$diff - mean(cov_2_table$diff)) / sd(cov_2_table$diff)
  fit_cov2 <- lm(diff_Zscaled ~ chr -1, data = cov_2_table)
  
  cov_3_table$diff_Zscaled <- (cov_3_table$diff - mean(cov_3_table$diff)) / sd(cov_3_table$diff)
  fit_cov3 <- lm(diff_Zscaled ~ chr -1, data = cov_3_table)
  
  snp_table$diff_Zscaled <- (snp_table$diff - mean(snp_table$diff)) / sd(snp_table$diff)
  fit_snp <- lm(diff_Zscaled ~ chr -1, data = snp_table)

  ################################################################################
  ############################# PLOTTING ################################
  ################################################################################

  text_size_colour = list(theme_bw(base_family="Courier", base_size = 12) + 
                            theme(axis.text.x= element_text(colour="black", size=12)) +
                            theme(axis.text.y= element_text(colour="black", size=12)) +
                            theme(axis.title.x = element_text(colour="black",size=15)) + 
                            theme(axis.title.y = element_text(colour="black",size=15)) + 
                            theme(plot.title=element_text(family="Courier", size=15, colour="black", hjust = 0.5) ) +
                            theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
                            theme(plot.margin= margin(1, 1, 1, 1, "mm")) +
                            theme(panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                                                  colour = "lightgrey"), 
                                  panel.grid.minor = element_line(size = 0.15, linetype = 'solid',
                                                                  colour = "lightgrey"))) 
  
  theme_title <- function(...) {
    theme_gray(base_family = "Courier") + 
      theme(plot.title = element_text(face = "bold"))
  }
  title_theme <- calc_element("plot.title", theme_title())
  
  theme_description <- function(...) {
    theme_gray(base_family = "Courier") + 
      theme(plot.title = element_text(face = "plain",size = 12, hjust = 0))
  }
  theme_description <- calc_element("plot.title", theme_description())

  plotLM <- plot_models(
    fit_cov1, fit_cov2, fit_cov3, fit_snp, vline.color = "black",
    m.labels = c(sprintf("Genome coverage %s", ED1), sprintf("Genome coverage %s", ED2), sprintf("Genome coverage %s", ED3), "Heterozygosity"),
    show.values = FALSE, show.p = FALSE, p.shape = TRUE, legend.title = "Measurement"
  ) + text_size_colour 
 
   data <- ggdraw() + draw_label(paste0("Data points from tables: \n", file1, " \n ",
                                        file2, " \n ", file3, " \n ", filesnp),     
                                 fontfamily = theme_description$family,
                                 fontface = theme_description$face,
                                 size = theme_description$size)
  
  data2 <- ggdraw() + draw_label(paste0("Results from linear models in HTML file: ", sprintf("%s.html", table_out_base)),     
    fontfamily = theme_description$family,
    fontface = theme_description$face,
    size = theme_description$size)
  
  title <- ggdraw() + draw_label("Linear Model Estimates and 95% CI",     
                                 fontfamily = title_theme$family,
                                 fontface = title_theme$face,
                                 size = title_theme$size)
  plotLM <- plot_grid(title, plotLM, nrow = 2, rel_heights = c(0.15, 1))

  data <- plot_grid(data, 
                    data2,
                   ncol = 1)

  len = length(unique(snp_table$chr))
  outname <- sprintf("%s.pdf", plot_out_base)
  pdf(file=outname, width = 18, height = len)
  print(plotLM)
  print(data)
  dev.off()  

  outname <- sprintf("%s.html", table_out_base)
  tab_model(fit_cov1, fit_cov2, fit_cov3, fit_snp, 
            title = "Results from linear models: Testing which chromosomes have significant differences between sexes (significantly different from 0)", 
            dv.labels = c(sprintf("Genome coverage %s", ED1), sprintf("Genome coverage %s", ED2), sprintf("Genome coverage %s", ED3), "Heterozygosity"),
            CSS = list(css.separatorcol = 'padding-right:2.5em; padding-left:2.5em;'), file = outname)
  
}
