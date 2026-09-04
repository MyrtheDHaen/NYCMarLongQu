################################################################################
# This file contains sample code for 'Quantile methods for longitudinal        #
# running data: a New York City Marathon example' (Under review, 2026) -       #
# D'Haen, Flórez, Van Keilegom, Molenberghs, Delecluse, Verhasselt.            #
#                                                                              #
# Code author: Myrthe D'Haen                                                   #
# Last revised on: 04/09/2026                                                  #
#                                                                              #
################################################################################



#### Part 0: initialisation (folders, functions, packages) #####################

# specify the correct working directory; the master directory should contain
# the data file NYClong.csv.
master.dir <- 'G:/My Drive/Marathon project'
setwd(master.dir)

library(rlist)
library(xtable)
library(ggplot2)
library(matrixStats)
library(patchwork)
library(RColorBrewer)
library(dplyr)
library(tidyverse)
library(skimr)
library(latex2exp)


# folders for output and plotting
store.EQCC.dir <- paste(master.dir, "EQCC numerical output", sep = "/")
plot.EQCC.dir <- paste(master.dir, "EQCC visual output", sep = "/")

plot.QB.dir <- paste(master.dir, "QB directory", sep = "/")
store.QB.EQCC.dir <- paste(plot.QB.dir, "EQCC numerical output", sep = "/")
plot.QB.EQCC.dir <- paste(plot.QB.dir, "EQCC visual output", sep = "/")
store.QB.EQCC.sens.dir <- paste(store.QB.EQCC.dir, "/Sensitivity", sep = "/")
plot.QB.EQCC.sens.dir <- paste(plot.QB.EQCC.dir, "/Sensitivity", sep = "/")


for (folder in c(store.EQCC.dir, plot.EQCC.dir, plot.QB.dir,
                 store.QB.EQCC.dir, plot.QB.EQCC.dir,
                 store.QB.EQCC.sens.dir, plot.QB.EQCC.sens.dir)){
  if (!dir.exists(folder)){
    dir.create(folder)
  }
}


## some basic functions --------------------------------------------------------
IntervalSpeed = function(data){ 
  c(data$dis[1]/(data$time[1]/60), diff(data$dis)/(diff(data$time)/60))
}
CumulativeSpeed = function(data){ 
  data$dis/(data$time/60)
}
SplittingVector = function(data){
  data$cumSpeed/tail(data$cumSpeed,1)
}
SpeedRatio = function(data){
  c(NA, tail(data$intSpeed, -1)/head(data$cumSpeed, -1))
}

## creating the "quantiles" of variable y.var (at each distance point)
# det.dist = OVERALL INDEX of the determining distance (for metrics that are
# distance-dependent, this determines which quantile will be used for logical
# subsetting). Indicates the position in the full distance vector (dist.vec),
# also containing any elements that may not be part of the selected dist.indices.

# NOTE: this function assumes that the supplied data contains only those rows
# corresponding to the distances provided (i.e., to those values in
# dist.vec[dist.indices]). Make sure any other rows are omitted before supplying
# the data frame to this function.
allQuantiles <- function(data, var.name, dist.var, dist.vec =
                           c(5,10,15,20,21.0975,25,30,32.1869,35,40,42.195),
                         dist.indices = 1:11, det.dist = NA){
  data$qu.lvl <- rep(0, dim(data)[1])
  data$qu.short <- rep(0, dim(data)[1])
  
  if (!is.na(det.dist)){
    if (!det.dist %in% dist.indices){
      print("Determining distance should be in the distance indices provided.")
      return(NA)
    }
  }
  
  for (index in seq_along(dist.indices)){
    this.data <- data[near(data[[dist.var]], dist.vec[dist.indices[index]]),]
    this.ecdf <- ecdf(x = this.data[[var.name]])
    
    this.qu.long <- this.ecdf(this.data[[var.name]])
    data[near(data[[dist.var]], dist.vec[dist.indices[index]]),]$qu.lvl <- this.qu.long
    
    # also append rounded version, that can be used for easier stratification
    # ensure that 0.020...01 up to 0.03 get assigned to 0.03 -> avoid 'round'
    this.qu.short <- 0.01*ceiling(100*this.qu.long)
    data[near(data[[dist.var]], dist.vec[dist.indices[index]]),]$qu.short <- this.qu.short
  }
  
  if (!is.na(det.dist)){
    data$qu.det <- rep(data$qu.short[near(data[[dist.var]], dist.vec[det.dist])],
                       each = length(dist.indices))
  } else { # default: final distance
    data$qu.det <- rep(data$qu.short[near(
      data[[dist.var]],tail(dist.vec[dist.indices],1))],
      each = length(dist.indices))
  }
  
  return(data)
}
# ------------------------------------------------------------------------------

################################################################################



#### Part 1: EQCC curves, with confidence intervals (CI) #######################

## Data preparation ------------------------------------------------------------
NYC.long = read_csv('NYClong.csv')
NYC.long <- mutate(NYC.long,
                   Sex = factor(Sex),
                   AgeGroup = factor(AgeGroup,
                                     levels = c("[18,19]","[20,29]","[30,39]",
                                                "[40,49]","[50,59]","[60,69]",
                                                "[70,79]", "[80,89]")),
                   distance = factor(distance))

NYC.long$Age2 = (NYC.long$Age - 18)/10
NYC.long$Mtime = rep(NYC.long$time[NYC.long$distance=="MAR"],each=11)

NYC.long$intSpeed = unlist(by(NYC.long,NYC.long$bib,IntervalSpeed))
NYC.long$cumSpeed = unlist(by(NYC.long,NYC.long$bib,CumulativeSpeed))
NYC.long$splitVec = unlist(by(NYC.long,NYC.long$bib,SplittingVector))
NYC.long$speedRat = unlist(by(NYC.long,NYC.long$bib,SpeedRatio))

NYC.longM <- NYC.long[NYC.long$Sex == "M",]
NYC.longF <- NYC.long[NYC.long$Sex == "F",]

## append quantiles
NYC.longM <- allQuantiles(NYC.longM, 'time', dist.var = "dis", dist.vec = 
                            c(5,10,15,20,21.0975,25,30,32.1869,35,40,42.195))
NYC.longF <- allQuantiles(NYC.longF, 'time', dist.var = "dis", dist.vec = 
                            c(5,10,15,20,21.0975,25,30,32.1869,35,40,42.195))

NYC.longM$ranking <- NYC.longM$qu.short
NYC.longF$ranking <- NYC.longF$qu.short
# ------------------------------------------------------------------------------


## Additional functions --------------------------------------------------------

## Auxiliary functions for logical indexing (manual implementation of equivalent
# to near(), of what %in% is to "==".)
near_in <- function(num, vec){
  for (index in seq_along(vec)){
    if (near(num, vec[index])){
      return(TRUE)
    }
  }
  return(FALSE)
}

# vectorized form
near_in_vec <- function(vec1,vec2){
  return(sapply(seq_along(vec1), function(ind1){
    near_in(num = vec1[ind1], vec = vec2)
  }))
}


## EQCC functions. Can also be applied to distance-dependent metrics. When
# applied to distance-independent metrics, just supply cond.dist = 11 or
# any other (present) distance; this choice is of no influence.

# function returning the curve values (without plotting, nor CI). Output in the
# form of a data frame with EQCC values per conditioning quantile level (rows)
# and distance index (columns).
get.EQCC <- function(df, y.var, qu.var, dist.var = "dis",
                     qu.lvls = c(0.25,0.5,0.75), dist.indices = (1:11),
                     dist.vec = c(5,10,15,20,21.0975,25,30,32.1869,35,40,42.195),
                     cond.dist = 11){
  # df = input data frame
  # y.var = name (string) of variable displayed on the EQCC curve
  # qu.var = name (string) of quantile-conditioning variable
  # dist.var = name (string) of variable containing the distances (in the same
  #       encoding as in dist.vec)
  # id.var = name (string) of id variable characterising individual runners
  # qu.lvls = levels at which quantile-conditioning is performed
  # dist.indices = measurement points that are to be included from dist.vec
  #       (often, either all points, or omitting only D1)
  # dist.vec = vector of all distances available in the supplied dataset (incl.
  #       any points where the EQCC values are not to be computed.)
  # cond.dist = index of the conditioning distance, at which the quantiles are
  #       determined and corresponding groups are selected. Should always be
  #       part of the dist.indices supplied.
  
  # extend data with all quantile levels wrt. the given variable, but select
  # only those distances (~rows) that are asked for.
  df.ext <- allQuantiles(data = df[
    near_in_vec(df[[dist.var]], dist.vec[dist.indices]),],
    var.name = qu.var, dist.var = dist.var, dist.vec = dist.vec,
    dist.indices = dist.indices, det.dist = cond.dist)
  
  # initialisation of output data frame
  df.EQCC <- data.frame(matrix(NA, nrow = length(qu.lvls),
                               ncol = length(dist.indices)))
  rownames(df.EQCC) <- paste0("qu.", qu.lvls)
  colnames(df.EQCC) <- paste0("D", dist.indices)
  
  for (qu.index in seq_along(qu.lvls)){
    this.qu = df.ext[((df.ext$qu.det > (qu.lvls[qu.index] - 0.001)) &
                        (df.ext$qu.det < (qu.lvls[qu.index] + 0.001))),]
    
    if (nrow(this.qu) > 0){
      # pointwise averages within this quantile
      for (index in seq_along(dist.indices)){
        # at each entry, store mean of all y.var values at that distance
        df.EQCC[qu.index, index] <- mean(this.qu[[y.var]][
          near(this.qu[[dist.var]], dist.vec[dist.indices[index]])])
      }
    } # else: preserve 'NA' from initialisation, since there are no observations
    # corresponding to this quantile level.
  }
  
  return(df.EQCC)
}


# wrapper function: compute EQCC curves and possibly CI bounds for one given df
EQCC.wrap <- function(df, y.var, qu.var, dist.var = "dis", id.var = "bib",
                      qu.lvls = c(0.25,0.5,0.75), dist.indices = (1:11),
                      dist.vec = c(5,10,15,20,21.0975,25,30,32.1869,35,40,42.195),
                      cond.dist = 11,
                      comp.CI = FALSE, store.path = NA,
                      sgnf = 0.05, boot.rep = 1000, bonf.corr = 1){
  
  # first arguments: see get.EQCC() documentation.
  
  # comp.CI = whether or not confidence interval is to be computed
  # store.path = folder where the EQCC values and any CI bounds are stored;
  #       put to NA if no external storage is desired.
  # sgnf = significance level for CI
  # boot.rep = number of replications used for bootstrapping
  # bonf.corr = the Bonferroni correction for multiple testing to be applied.
  #       Default value = 1 corresponds to no correction; in general a p-value
  #       of sgnf/bonf.corr is used.
  
  # Note that CI are no full confidence bands, but rather pointwise CI computed
  # at each distance included in dist.indices.
  
  ### Computation of EQCC values
  df.EQCC <- get.EQCC(df = df, y.var = y.var, qu.var = qu.var,
                      dist.var = dist.var, qu.lvls = qu.lvls,
                      dist.indices = dist.indices, dist.vec = dist.vec,
                      cond.dist = cond.dist)
  
  ### Computation of EQCC confidence intervals
  if (comp.CI){
    
    ## Compute EQCC values for boot.rep bootstrap resamples of the data
    
    # Initialisation: list of data frames for storage of EQCC values for each
    # bootstrap resample; rows correspond to quantile levels, cols to distances.
    df.empty <- data.frame(matrix(NA, nrow = length(qu.lvls),
                                  ncol = length(dist.indices)))
    rownames(df.empty) <- paste0("qu.", qu.lvls)
    colnames(df.empty) <- paste0("D", dist.indices)
    df.list <- replicate(boot.rep, df.empty, simplify = FALSE)
    
    all.ids <- unique(df[[id.var]])
    
    # for Bonferroni correction
    sgnf <- sgnf/bonf.corr
    
    for (boot.index in seq_len(boot.rep)){
      
      ## printing progress
      if (boot.index %% 25 == 0){
        print(boot.index)
      }
      
      ## bootstrap resample
      set.seed(boot.index)
      resample.ids <- sample(all.ids, size = length(all.ids), replace = TRUE)
      
      # avoid logical indexing, where items selected multiple times appear only once
      df.bibs <- data.frame(matrix(resample.ids, ncol = 1))
      colnames(df.bibs) <- id.var
      
      if (id.var == "bib"){ # avoid 'join_by' message by explicitly including it
        df.resampled <- inner_join(x = df, y = df.bibs,
                                   relationship = "many-to-many", by = join_by(bib))
      } else { # also works for id.var = "bib", but throws a message each time
        df.resampled <- inner_join(x = df, y = df.bibs,
                                   relationship = "many-to-many")
      }
      
      
      ## EQCC computation
      df.list[[boot.index]] <- get.EQCC(
        df = df.resampled, y.var = y.var, qu.var = qu.var, dist.var = dist.var,
        qu.lvls = qu.lvls, dist.indices = dist.indices, dist.vec = dist.vec,
        cond.dist = cond.dist)
    }
    
    ## Compute lower and upper bounds corresponding to the given significance
    
    # Initialisation
    df.lower <- df.empty
    df.upper <- df.empty
    
    for (qu.index in seq_along(qu.lvls)){
      for (dist.index in seq_along(dist.indices)){
        # combine all (boot.rep) computed values at this entry, into a vector
        vals <- sapply(seq_len(boot.rep), function(boot.index){
          return(df.list[[boot.index]][qu.index,dist.index])
        })
        
        
        # take care of NA values that may occur (as some quantile levels may
        # no longer be present in certain bootstrap resamples)
        if (any(is.na(vals))){
          vals.noNA <- vals[!is.na(vals)]
          vals.sorted <- sort(vals.noNA)
          rep.remaining <- length(vals.noNA)
          
          # Throw a warning only once; in case any NA values occur, they always
          # occur at each distance (cf. get.EQCC function: entire row of NA then)
          # (but different per qu.lvl -> warning for each qu.index where relevant)
          # (For monitoring that not too many resamples are omitted.)
          if (dist.index == 1){
            if (rep.remaining == (boot.rep - 1)){
              warning(paste0("At level qu = ", qu.lvls[qu.index],
                             ", 1 out of ", boot.rep, " bootstrap resamples",
                             " was omitted due to NA values."))
            } else {
              warning(paste0("At level qu = ", qu.lvls[qu.index], ", ",
                             boot.rep - rep.remaining,
                             " out of ", boot.rep, " bootstrap resamples",
                             " were omitted due to NA values."))
            }
          }
        } else {
          vals.sorted <- sort(vals)
          rep.remaining <- boot.rep
        }
        
        
        # select lower and upper bounds such that only proportion (at most)
        # sgnf/2 below and above these values, respectively.
        df.lower[qu.index, dist.index] <-
          vals.sorted[max(floor(rep.remaining*sgnf/2),1)]
        # max() as safety catch for too small boot.rep, where floor(...) = 0
        df.upper[qu.index, dist.index] <-
          vals.sorted[ceiling(rep.remaining*(1-sgnf/2))]
      }
    }
    
  } else { # not actual data frames, but to avoid case distinction below
    df.lower <- NA
    df.upper <- NA
  }
  
  ### Storage of results
  if (!is.na(store.path)){
    
    if (! dir.exists(store.path)){
      dir.create(store.path)
    }
    
    save(df.EQCC, file = paste0(store.path, "/df.", y.var ,".EQCC.Rdata"))
    
    if (comp.CI){
      save(df.lower, file = paste0(store.path, "/df.", y.var ,".CI.lower.Rdata"))
      save(df.upper, file = paste0(store.path, "/df.", y.var ,".CI.upper.Rdata"))
    }
  }
  
  
  ### Also return results to console
  return(list(EQCC = df.EQCC, lower = df.lower, upper = df.upper))
}

# function for plotting EQCC. Possibly including CI, possibly combining multiple
# lists (male, female, different performance groups) to get one common scale.
EQCC.plot <- function(EQCC.out.list, y.var, qu.var, y.label, title.vec,
                      storeName.vec, box.labels, box.locations,
                      dist.var = "dis", id.var = "bib",
                      qu.lvls = c(0.25,0.5,0.75), dist.indices = (1:11),
                      dist.vec = c(5,10,15,20,21.0975,25,30,32.1869,35,40,42.195),
                      plot.CI = FALSE, store.path = NA,
                      save.separate.panels = FALSE){
  
  # basic arguments: see get.EQCC() and EQCC.wrap() documentation.
  
  # EQCC.out.list = list containing EQCC output lists (as produced by EQCC.wrap),
  #       that are to be visualised. Each entry should be a named list with sub-
  #       entries $EQCC, $lower, and $upper.
  # y.label = string to be used as axis label for the y.var variable
  # title.vec = vector (length as df.list) of strings to be used as plot titles
  # storeName.vec = vector (same length as df.list) of strings to be used as
  #       name strings to store the separate panels
  # NOTE: title.vec and storeName.vec can also contain one additional element at
  #       the end, that contains the strings to be used for the overview plot
  #       combining all panels. If not provided, the first elements are re-used.
  # box.labels = vector (same length as df.list) containing labels to be put in
  #       the 'legend' boxes (typically of the following form, possibly
  #       supplemented with qu.var: c(Female, Female, Female, Male, Male, Male))
  # box.locations = list (same length as df.list) containing coordinate tuples
  #       indicating where to put the 'legend' boxes.
  # plot.CI = whether or not a CI is to be plotted
  # save.separate.panels = whether or not all panels are to be stored as separate
  #       plots; they are always arranged and stored in one overall plot as well.
  
  
  
  ## prepare dfs for more convenient ggplotting
  if (plot.CI){ # case with CI
    
    # verifying whether everything is properly supplied; mostly serves to catch
    # the problem of putting plot.CI to TRUE while no CI bounds are supplied.
    if (anyNA(EQCC.out.list, recursive = TRUE)){
      print(paste("The supplied arguments contain NA values. Please make sure",
                  "that any lower and upper bounds to be plotted, are supplied."))
      return(NA)
    }
    
    dfs.prepared <- lapply(seq_along(EQCC.out.list), function(df.index){
      df.prep <- data.frame(cbind(rep(qu.lvls, each =
                                        length(dist.indices)),
                                  rep(dist.vec[dist.indices],length(qu.lvls)),
                                  as.vector(t(EQCC.out.list[[df.index]]$EQCC)),
                                  as.vector(t(EQCC.out.list[[df.index]]$lower)),
                                  as.vector(t(EQCC.out.list[[df.index]]$upper))))
      colnames(df.prep) <- c("qu.lvl", "distance", "EQCC.val", "EQCC.lower",
                             "EQCC.upper")
      return(df.prep)
    })
    
    # common range over all data subsets, to ensure plotting on the same scale
    y.range <- range(EQCC.out.list)
  } else { # case without CI
    dfs.prepared <- lapply(seq_along(EQCC.out.list), function(df.index){
      df.prep <- data.frame(cbind(rep(qu.lvls, each =
                                        length(dist.indices)),
                                  rep(dist.vec[dist.indices],length(qu.lvls)),
                                  as.vector(t(EQCC.out.list[[df.index]]$EQCC))))
      colnames(df.prep) <- c("qu.lvl", "distance", "EQCC.val")
      return(df.prep)
    })
    
    # common range over all data subsets, to ensure plotting on the same scale
    y.range <- range(EQCC.out.list)
  }
  
  ## preparing plot panels
  all.panels <- vector(mode = "list", length = length(EQCC.out.list))
  for (plot.index in seq_along(EQCC.out.list)){
    df.prep <- dfs.prepared[[plot.index]]
    
    if (plot.CI){ # case with CI
      p <- ggplot(df.prep, aes(x = distance, y = EQCC.val,
                               group = factor(qu.lvl), color = factor(qu.lvl))) +
        geom_line(linewidth = 1.2) + geom_point(size = 1.7) +
        geom_ribbon(aes(ymin = EQCC.lower, ymax = EQCC.upper,
                        fill = factor(qu.lvl), color = factor(qu.lvl)),
                    alpha = 0.3) +
        xlab('Distance (km)') + ylab(y.label) + ylim(y.range) +
        scale_colour_manual(name = "Levels", values = brewer.pal(length(qu.lvls),
                                                                 "YlGnBu")) +
        scale_fill_manual(name = "Levels", values = brewer.pal(length(qu.lvls),
                                                               "YlGnBu")) +
        theme(legend.position = "right") + theme(text = element_text(size = 15)) +
        theme(panel.background = element_rect(fill = "grey90", colour = "grey90",
                                              linewidth = 2, linetype = "solid")) +
        geom_label(label = box.labels[plot.index],
                   x = box.locations[[plot.index]][1],
                   y = box.locations[[plot.index]][2],
                   label.padding = unit(0.40, "lines"),
                   linewidth = 0.35, colour = "black", fill = "white")
    } else { # case without CI
      p <- ggplot(df.prep, aes(x = distance, y = EQCC.val,
                               group = factor(qu.lvl), color = factor(qu.lvl))) +
        geom_line(linewidth = 1.2) + geom_point(size = 1.7) +
        xlab('Distance (km)') + ylab(y.label) + ylim(y.range) +
        scale_colour_manual(name = "Levels", values = brewer.pal(length(qu.lvls),
                                                                 "YlGnBu")) +
        theme(legend.position = "right") + theme(text = element_text(size = 15)) +
        theme(panel.background = element_rect(fill = "grey90", colour = "grey90",
                                              linewidth = 2, linetype = "solid")) +
        geom_label(label = box.labels[plot.index],
                   x = box.locations[[plot.index]][1],
                   y = box.locations[[plot.index]][2],
                   label.padding = unit(0.40, "lines"),
                   linewidth = 0.35, colour = "black", fill = "white")
    }
    
    all.panels[[plot.index]] <- p
    
    if (save.separate.panels){
      p.titled <- p + ggtitle(title.vec[plot.index])
      ggsave(path = store.path,
             filename = paste0(storeName.vec[plot.index], ".png"),
             width = 6, height = 4, device = "png")
    }
  }
  
  ## Combining plot panels
  # preparation of title and file name
  if (length(title.vec) > length(EQCC.out.list)){
    title.combi <- title.vec[length(EQCC.out.list) + 1]
  } else {
    title.combi <- title.vec[1]
  }
  
  if (length(storeName.vec) > length(EQCC.out.list)){
    storeName.combi <- storeName.vec[length(EQCC.out.list) + 1]
  } else {
    storeName.combi <- paste0(storeName.vec[1], "_Combi")
  }
  
  # arranging all panels (note: by column by default, but can be modified here)
  pCombi <- wrap_plots(all.panels, ncol = 2, byrow = FALSE) +
    plot_layout(guides = 'collect', axes = 'collect_x',
                axis_titles = 'collect') +
    plot_annotation(title = title.combi,
                    theme = theme(plot.title = element_text(hjust = 0.47),
                                  text = element_text(size=15)))
  ggsave(path = store.path,
         filename = paste0(storeName.combi,".png"),
         width = 10, height = 4*ceiling(length(all.panels)/2), device='png')
}


## faster alternative to EQCC.wrap (as this takes long, especially with 1000
# bootstrap replications), but tailored to the specific analysis performed. The
# idea of incorporating all variables of interest into one function (since
# bootstrap resampling and extending with quantiles based on "time" are common
# over all possible y.var) can, however, be generalised.
# One can also incorporate both a Bonferroni corrected and non-corrected version.
# It is moreover convenient to specify as many relevant quantile levels (qu.lvls)
# as possible to avoid having to rerun all code for addition of an extra level.
# Unused quantile levels can simply be omitted for the plotting later.
EQCC.multi <- function(df, qu.var, dist.var = "dis", id.var = "bib",
                       qu.lvls = c(0.01,0.05,0.1,0.2,0.25,0.3,0.4,0.5,
                                   0.6,0.7,0.75,0.8,0.9,0.95,0.99),
                       dist.vec = c(5,10,15,20,21.0975,25,30,32.1869,35,40,42.195),
                       cond.dist = 11,
                       comp.CI = FALSE, store.path = NA, store.name = NA,
                       sgnf = 0.05, boot.rep = 1000){
  
  ### Computation of EQCC values
  df.EQCC.intSpeed <- get.EQCC(df = df, y.var = "intSpeed", qu.var = qu.var,
                               dist.var = dist.var, qu.lvls = qu.lvls,
                               dist.indices = 1:11, dist.vec = dist.vec,
                               cond.dist = cond.dist)
  df.EQCC.cumSpeed <- get.EQCC(df = df, y.var = "cumSpeed", qu.var = qu.var,
                               dist.var = dist.var, qu.lvls = qu.lvls,
                               dist.indices = 1:11, dist.vec = dist.vec,
                               cond.dist = cond.dist)
  df.EQCC.splitVec <- get.EQCC(df = df, y.var = "splitVec", qu.var = qu.var,
                               dist.var = dist.var, qu.lvls = qu.lvls,
                               dist.indices = 1:11, dist.vec = dist.vec,
                               cond.dist = cond.dist)
  df.EQCC.speedRat <- get.EQCC(df = df, y.var = "speedRat", qu.var = qu.var,
                               dist.var = dist.var, qu.lvls = qu.lvls,
                               dist.indices = 2:11, dist.vec = dist.vec,
                               cond.dist = cond.dist)
  
  
  ### Computation of EQCC confidence intervals
  if (comp.CI){
    
    ## Compute EQCC values for boot.rep bootstrap resamples of the data
    
    # Initialisation: list of data frames for storage of EQCC values for each
    # bootstrap resample; rows correspond to quantile levels, cols to distances.
    df.empty.11 <- data.frame(matrix(NA, nrow = length(qu.lvls), ncol = 11))
    rownames(df.empty.11) <- paste0("qu.", qu.lvls)
    colnames(df.empty.11) <- paste0("D", 1:11)
    
    df.empty.10 <- data.frame(matrix(NA, nrow = length(qu.lvls), ncol = 10))
    rownames(df.empty.10) <- paste0("qu.", qu.lvls)
    colnames(df.empty.10) <- paste0("D", 2:11)
    
    df.list.intSpeed <- replicate(boot.rep, df.empty.11, simplify = FALSE)
    df.list.cumSpeed <- replicate(boot.rep, df.empty.11, simplify = FALSE)
    df.list.splitVec <- replicate(boot.rep, df.empty.11, simplify = FALSE)
    df.list.speedRat <- replicate(boot.rep, df.empty.10, simplify = FALSE)
    
    all.ids <- unique(df[[id.var]])
    
    for (boot.index in seq_len(boot.rep)){
      
      ## printing progress
      if (boot.index %% 25 == 0){
        print(boot.index)
      }
      
      ## bootstrap resample
      set.seed(boot.index)
      resample.ids <- sample(all.ids, size = length(all.ids), replace = TRUE)
      
      # avoid logical indexing, where items selected multiple times appear only once
      df.bibs <- data.frame(matrix(resample.ids, ncol = 1))
      colnames(df.bibs) <- id.var
      
      if (id.var == "bib"){ # avoid 'join_by' message by explicitly including it
        df.resampled <- inner_join(x = df, y = df.bibs,
                                   relationship = "many-to-many", by = join_by(bib))
      } else { # also works for id.var = "bib", but throws a message each time
        df.resampled <- inner_join(x = df, y = df.bibs,
                                   relationship = "many-to-many")
      }
      
      ## EQCC computation
      df.list.intSpeed[[boot.index]] <- get.EQCC(
        df = df.resampled, y.var = "intSpeed", qu.var = qu.var, dist.var = dist.var,
        qu.lvls = qu.lvls, dist.indices = 1:11, dist.vec = dist.vec,
        cond.dist = cond.dist)
      df.list.cumSpeed[[boot.index]] <- get.EQCC(
        df = df.resampled, y.var = "cumSpeed", qu.var = qu.var, dist.var = dist.var,
        qu.lvls = qu.lvls, dist.indices = 1:11, dist.vec = dist.vec,
        cond.dist = cond.dist)
      df.list.splitVec[[boot.index]] <- get.EQCC(
        df = df.resampled, y.var = "splitVec", qu.var = qu.var, dist.var = dist.var,
        qu.lvls = qu.lvls, dist.indices = 1:11, dist.vec = dist.vec,
        cond.dist = cond.dist)
      df.list.speedRat[[boot.index]] <- get.EQCC(
        df = df.resampled, y.var = "speedRat", qu.var = qu.var, dist.var = dist.var,
        qu.lvls = qu.lvls, dist.indices = 2:11, dist.vec = dist.vec,
        cond.dist = cond.dist)
    }
    
    ## Compute lower and upper bounds corresponding to the given significance
    
    # Initialisation
    df.lower.intSpeed <- df.empty.11
    df.upper.intSpeed <- df.empty.11
    
    df.lower.cumSpeed <- df.empty.11
    df.upper.cumSpeed <- df.empty.11
    
    df.lower.splitVec <- df.empty.11
    df.upper.splitVec <- df.empty.11
    
    df.lower.speedRat <- df.empty.10
    df.upper.speedRat <- df.empty.10
    
    for (qu.index in seq_along(qu.lvls)){
      for (dist.index in (1:11)){
        
        vals.intSpeed <- sapply(seq_len(boot.rep), function(boot.index){
          return(df.list.intSpeed[[boot.index]][qu.index,dist.index])
        })
        vals.cumSpeed <- sapply(seq_len(boot.rep), function(boot.index){
          return(df.list.cumSpeed[[boot.index]][qu.index,dist.index])
        })
        vals.splitVec <- sapply(seq_len(boot.rep), function(boot.index){
          return(df.list.splitVec[[boot.index]][qu.index,dist.index])
        })
        
        
        # take care of any NA values (note: rows with NA values are identical
        # for all variables considered, as they are determined by qu.var, that
        # is common over all y.var)
        if (any(is.na(vals.intSpeed))){
          
          vals.intSpeed.noNA <- vals.intSpeed[!is.na(vals.intSpeed)]
          vals.intSpeed.sorted <- sort(vals.intSpeed.noNA)
          
          vals.cumSpeed.noNA <- vals.cumSpeed[!is.na(vals.cumSpeed)]
          vals.cumSpeed.sorted <- sort(vals.cumSpeed.noNA)
          
          vals.splitVec.noNA <- vals.splitVec[!is.na(vals.splitVec)]
          vals.splitVec.sorted <- sort(vals.splitVec.noNA)
          
          rep.remaining <- length(vals.intSpeed.noNA)
          
          if (dist.index == 1){
            if (rep.remaining == (boot.rep - 1)){
              warning(paste0("At level qu = ", qu.lvls[qu.index],
                             ", 1 out of ", boot.rep, " bootstrap resamples",
                             " was omitted due to NA values."))
            } else {
              warning(paste0("At level qu = ", qu.lvls[qu.index], ", ",
                             boot.rep - rep.remaining,
                             " out of ", boot.rep, " bootstrap resamples",
                             " were omitted due to NA values."))
            }
          }
        } else {
          vals.intSpeed.sorted <- sort(vals.intSpeed)
          vals.cumSpeed.sorted <- sort(vals.cumSpeed)
          vals.splitVec.sorted <- sort(vals.splitVec)
          rep.remaining <- boot.rep
        }
        
        
        
        df.lower.intSpeed[qu.index, dist.index] <-
          vals.intSpeed.sorted[max(floor(rep.remaining*sgnf/2),1)]
        df.upper.intSpeed[qu.index, dist.index] <-
          vals.intSpeed.sorted[ceiling(rep.remaining*(1-sgnf/2))]
        
        df.lower.cumSpeed[qu.index, dist.index] <-
          vals.cumSpeed.sorted[max(floor(rep.remaining*sgnf/2),1)]
        df.upper.cumSpeed[qu.index, dist.index] <-
          vals.cumSpeed.sorted[ceiling(rep.remaining*(1-sgnf/2))]
        
        df.lower.splitVec[qu.index, dist.index] <-
          vals.splitVec.sorted[max(floor(rep.remaining*sgnf/2),1)]
        df.upper.splitVec[qu.index, dist.index] <-
          vals.splitVec.sorted[ceiling(rep.remaining*(1-sgnf/2))]
        
        
        if (dist.index > 1){ # speedRat defined only for dist.indices = 2:11
          vals.speedRat <- sapply(seq_len(boot.rep), function(boot.index){
            return(df.list.speedRat[[boot.index]][qu.index,dist.index-1])
          })
          
          if (any(is.na(vals.speedRat))){
            vals.speedRat.noNA <- vals.speedRat[!is.na(vals.speedRat)]
            vals.speedRat.sorted <- sort(vals.speedRat.noNA)
          } else {
            vals.speedRat.sorted <- sort(vals.speedRat)
          }
          
          df.lower.speedRat[qu.index, dist.index-1] <-
            vals.speedRat.sorted[max(floor(rep.remaining*sgnf/2),1)]
          df.upper.speedRat[qu.index, dist.index-1] <-
            vals.speedRat.sorted[ceiling(rep.remaining*(1-sgnf/2))]
        }
      }
    }
    
  } else { # not actual data frames, but to avoid case distinction below
    df.lower.intSpeed <- NA
    df.upper.intSpeed <- NA
    
    df.lower.cumSpeed <- NA
    df.upper.cumSpeed <- NA
  
    df.lower.splitVec <- NA
    df.upper.splitVec <- NA
    
    df.lower.speedRat <- NA
    df.upper.speedRat <- NA
  }
  
  
  full.output <- list(
    resIntSpeed = list(EQCC = df.EQCC.intSpeed,
                       lower = df.lower.intSpeed, upper = df.upper.intSpeed),
    resCumSpeed = list(EQCC = df.EQCC.cumSpeed,
                       lower = df.lower.cumSpeed, upper = df.upper.cumSpeed),
    resSplitVec = list(EQCC = df.EQCC.splitVec,
                       lower = df.lower.splitVec, upper = df.upper.splitVec),
    resSpeedRat = list(EQCC = df.EQCC.speedRat,
                       lower = df.lower.speedRat, upper = df.upper.speedRat)
  )
  
  ### Storage of results
  if (!is.na(store.path)){
    
    if (! dir.exists(store.path)){
      dir.create(store.path)
    }
    
    list.save(full.output,
              file = paste0(store.path, "/full.output.", store.name, ".Rdata"))
  }
  
  
  ### Also return results to console
  return(full.output)
}
# ------------------------------------------------------------------------------


## Application of these functions ----------------------------------------------
# (alternatively, apply EQCC.wrap() to the variables of interest separately)

# Computation & (external) storage of the results, with conditioning distances
# d2 = 10 km, d5 = half marathon, d11 = full marathon

EQCC.M.d2.multi <- EQCC.multi(df = NYC.longM, qu.var = "time",
                              cond.dist = 2, comp.CI = TRUE,
                              store.path = store.EQCC.dir,
                              store.name = "M-d2",
                              boot.rep = 1000)
EQCC.M.d5.multi <- EQCC.multi(df = NYC.longM, qu.var = "time",
                              cond.dist = 5, comp.CI = TRUE,
                              store.path = store.EQCC.dir,
                              store.name = "M-d5",
                              boot.rep = 1000)
EQCC.M.d11.multi <- EQCC.multi(df = NYC.longM, qu.var = "time",
                               cond.dist = 11, comp.CI = TRUE,
                               store.path = store.EQCC.dir,
                               store.name = "M-d11",
                               boot.rep = 1000)

EQCC.F.d2.multi <- EQCC.multi(df = NYC.longF, qu.var = "time",
                              cond.dist = 2, comp.CI = TRUE,
                              store.path = store.EQCC.dir,
                              store.name = "F-d2",
                              boot.rep = 1000)
EQCC.F.d5.multi <- EQCC.multi(df = NYC.longF, qu.var = "time",
                              cond.dist = 5, comp.CI = TRUE,
                              store.path = store.EQCC.dir,
                              store.name = "F-d5",
                              boot.rep = 1000)
EQCC.F.d11.multi <- EQCC.multi(df = NYC.longF, qu.var = "time",
                               cond.dist = 11, comp.CI = TRUE,
                               store.path = store.EQCC.dir,
                               store.name = "F-d11",
                               boot.rep = 1000)
# ------------------------------------------------------------------------------


### Conversion to proper format, and plotting ----------------------------------

## Select quantile levels to be displayed
all.qu.lvls <- # default for EQCC.multi
  c(0.01,0.05,0.1,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.75,0.8,0.9,0.95,0.99)
qu.sel <- c(1,3,8,13,15) # for (0.01,0.1,0.5,0.9,0.99); other options possible

## retrieving data frames per variable ----
# (can be skipped if EQCC.wrap was applied per variable, rather than EQCC.multi)

# # in case they were computed in a previous session
# EQCC.M.d2.multi <- list.load(paste0(store.EQCC.dir,"/full.output.M-d2.Rdata"))
# EQCC.M.d5.multi <- list.load(paste0(store.EQCC.dir,"/full.output.M-d5.Rdata"))
# EQCC.M.d11.multi <- list.load(paste0(store.EQCC.dir,"/full.output.M-d11.Rdata"))
# EQCC.F.d2.multi <- list.load(paste0(store.EQCC.dir,"/full.output.F-d2.Rdata"))
# EQCC.F.d5.multi <- list.load(paste0(store.EQCC.dir,"/full.output.F-d5.Rdata"))
# EQCC.F.d11.multi <- list.load(paste0(store.EQCC.dir,"/full.output.F-d11.Rdata"))

intSpeed.M.d2 <- EQCC.M.d2.multi$resIntSpeed
intSpeed.M.d5 <- EQCC.M.d5.multi$resIntSpeed
intSpeed.M.d11 <- EQCC.M.d11.multi$resIntSpeed
intSpeed.F.d2 <- EQCC.F.d2.multi$resIntSpeed
intSpeed.F.d5 <- EQCC.F.d5.multi$resIntSpeed
intSpeed.F.d11 <- EQCC.F.d11.multi$resIntSpeed

cumSpeed.M.d2 <- EQCC.M.d2.multi$resCumSpeed
cumSpeed.M.d5 <- EQCC.M.d5.multi$resCumSpeed
cumSpeed.M.d11 <- EQCC.M.d11.multi$resCumSpeed
cumSpeed.F.d2 <- EQCC.F.d2.multi$resCumSpeed
cumSpeed.F.d5 <- EQCC.F.d5.multi$resCumSpeed
cumSpeed.F.d11 <- EQCC.F.d11.multi$resCumSpeed

splitVec.M.d2 <- EQCC.M.d2.multi$resSplitVec
splitVec.M.d5 <- EQCC.M.d5.multi$resSplitVec
splitVec.M.d11 <- EQCC.M.d11.multi$resSplitVec
splitVec.F.d2 <- EQCC.F.d2.multi$resSplitVec
splitVec.F.d5 <- EQCC.F.d5.multi$resSplitVec
splitVec.F.d11 <- EQCC.F.d11.multi$resSplitVec

speedRat.M.d2 <- EQCC.M.d2.multi$resSpeedRat
speedRat.M.d5 <- EQCC.M.d5.multi$resSpeedRat
speedRat.M.d11 <- EQCC.M.d11.multi$resSpeedRat
speedRat.F.d2 <- EQCC.F.d2.multi$resSpeedRat
speedRat.F.d5 <- EQCC.F.d5.multi$resSpeedRat
speedRat.F.d11 <- EQCC.F.d11.multi$resSpeedRat
# ----

## interval speed ----
v.M2 <- list(EQCC = intSpeed.M.d2$EQCC[qu.sel,],
             lower = intSpeed.M.d2$lower[qu.sel,],
             upper = intSpeed.M.d2$upper[qu.sel,])
v.M5 <- list(EQCC = intSpeed.M.d5$EQCC[qu.sel,],
             lower = intSpeed.M.d5$lower[qu.sel,],
             upper = intSpeed.M.d5$upper[qu.sel,])
v.M11 <- list(EQCC = intSpeed.M.d11$EQCC[qu.sel,],
              lower = intSpeed.M.d11$lower[qu.sel,],
              upper = intSpeed.M.d11$upper[qu.sel,])

v.F2 <- list(EQCC = intSpeed.F.d2$EQCC[qu.sel,],
             lower = intSpeed.F.d2$lower[qu.sel,],
             upper = intSpeed.F.d2$upper[qu.sel,])
v.F5 <- list(EQCC = intSpeed.F.d5$EQCC[qu.sel,],
             lower = intSpeed.F.d5$lower[qu.sel,],
             upper = intSpeed.F.d5$upper[qu.sel,])
v.F11 <- list(EQCC = intSpeed.F.d11$EQCC[qu.sel,],
              lower = intSpeed.F.d11$lower[qu.sel,],
              upper = intSpeed.F.d11$upper[qu.sel,])

EQCC.plot(EQCC.out.list = list(v.F2,v.F5,v.F11,v.M2,v.M5,v.M11),
          y.var = "intSpeed", qu.var = "time",
          y.label = "Interval speed (km/h)",
          title.vec = c("EQCC interval speed (F: 10K)", "EQCC interval speed (F: Half)",
                        "EQCC interval speed (F: Full)", "EQCC interval speed (M: 10K)",
                        "EQCC interval speed (M: Half)", "EQCC interval speed (M: Full)",
                        "EQCC interval speed"),
          storeName.vec = c("intSpeed-F-d2", "intSpeed-F-d5", "intSpeed-F-d11",
                            "intSpeed-M-d2", "intSpeed-M-d5", "intSpeed-M-d11",
                            "intSpeed-all"),
          box.labels = c("F: 10K", "F: Half", "F: Full",
                         "M: 10K", "M: Half", "M: Full"),
          box.locations = replicate(6,c(40,17.1), simplify = FALSE),
          qu.lvls = all.qu.lvls[qu.sel],
          dist.indices = 1:11,
          store.path = plot.EQCC.dir,
          plot.CI = TRUE, save.separate.panels = FALSE)
# ----

## cumulative speed ----
c.M2 <- list(EQCC = cumSpeed.M.d2$EQCC[qu.sel,],
             lower = cumSpeed.M.d2$lower[qu.sel,],
             upper = cumSpeed.M.d2$upper[qu.sel,])
c.M5 <- list(EQCC = cumSpeed.M.d5$EQCC[qu.sel,],
             lower = cumSpeed.M.d5$lower[qu.sel,],
             upper = cumSpeed.M.d5$upper[qu.sel,])
c.M11 <- list(EQCC = cumSpeed.M.d11$EQCC[qu.sel,],
              lower = cumSpeed.M.d11$lower[qu.sel,],
              upper = cumSpeed.M.d11$upper[qu.sel,])

c.F2 <- list(EQCC = cumSpeed.F.d2$EQCC[qu.sel,],
             lower = cumSpeed.F.d2$lower[qu.sel,],
             upper = cumSpeed.F.d2$upper[qu.sel,])
c.F5 <- list(EQCC = cumSpeed.F.d5$EQCC[qu.sel,],
             lower = cumSpeed.F.d5$lower[qu.sel,],
             upper = cumSpeed.F.d5$upper[qu.sel,])
c.F11 <- list(EQCC = cumSpeed.F.d11$EQCC[qu.sel,],
              lower = cumSpeed.F.d11$lower[qu.sel,],
              upper = cumSpeed.F.d11$upper[qu.sel,])

EQCC.plot(EQCC.out.list = list(c.F2,c.F5,c.F11,c.M2,c.M5,c.M11),
          y.var = "cumSpeed", qu.var = "time",
          y.label = "Cumulative speed (km/h)",
          title.vec = c("EQCC cumulative speed (F: 10K)", "EQCC cumulative speed (F: Half)",
                        "EQCC cumulative speed (F: Full)", "EQCC cumulative speed (M: 10K)",
                        "EQCC cumulative speed (M: Half)", "EQCC cumulative speed (M: Full)",
                        "EQCC cumulative speed"),
          storeName.vec = c("cumSpeed-F-d2", "cumSpeed-F-d5", "cumSpeed-F-d11",
                            "cumSpeed-M-d2", "cumSpeed-M-d5", "cumSpeed-M-d11",
                            "cumSpeed-all"),
          box.labels = c("F: 10K", "F: Half", "F: Full",
                         "M: 10K", "M: Half", "M: Full"),
          box.locations = replicate(6,c(40,17.1), simplify = FALSE),
          qu.lvls = all.qu.lvls[qu.sel],
          dist.indices = 1:11,
          store.path = plot.EQCC.dir,
          plot.CI = TRUE, save.separate.panels = FALSE)
# ----

## splitting indicator ----
s.M2 <- list(EQCC = splitVec.M.d2$EQCC[qu.sel,],
             lower = splitVec.M.d2$lower[qu.sel,],
             upper = splitVec.M.d2$upper[qu.sel,])
s.M5 <- list(EQCC = splitVec.M.d5$EQCC[qu.sel,],
             lower = splitVec.M.d5$lower[qu.sel,],
             upper = splitVec.M.d5$upper[qu.sel,])
s.M11 <- list(EQCC = splitVec.M.d11$EQCC[qu.sel,],
              lower = splitVec.M.d11$lower[qu.sel,],
              upper = splitVec.M.d11$upper[qu.sel,])

s.F2 <- list(EQCC = splitVec.F.d2$EQCC[qu.sel,],
             lower = splitVec.F.d2$lower[qu.sel,],
             upper = splitVec.F.d2$upper[qu.sel,])
s.F5 <- list(EQCC = splitVec.F.d5$EQCC[qu.sel,],
             lower = splitVec.F.d5$lower[qu.sel,],
             upper = splitVec.F.d5$upper[qu.sel,])
s.F11 <- list(EQCC = splitVec.F.d11$EQCC[qu.sel,],
              lower = splitVec.F.d11$lower[qu.sel,],
              upper = splitVec.F.d11$upper[qu.sel,])

EQCC.plot(EQCC.out.list = list(s.F2,s.F5,s.F11,s.M2,s.M5,s.M11),
          y.var = "splitVec", qu.var = "time",
          y.label = "Splitting indicator",
          title.vec = c("EQCC splitting indicator (F: 10K)", "EQCC splitting indicator (F: Half)",
                        "EQCC splitting indicator (F: Full)", "EQCC splitting indicator (M: 10K)",
                        "EQCC splitting indicator (M: Half)", "EQCC splitting indicator (M: Full)",
                        "EQCC splitting indicator"),
          storeName.vec = c("splitVec-F-d2", "splitVec-F-d5", "splitVec-F-d11",
                            "splitVec-M-d2", "splitVec-M-d5", "splitVec-M-d11",
                            "splitVec-all"),
          box.labels = c("F: 10K", "F: Half", "F: Full",
                         "M: 10K", "M: Half", "M: Full"),
          box.locations = replicate(6,c(40,1.275), simplify = FALSE),
          qu.lvls = all.qu.lvls[qu.sel],
          dist.indices = 1:11,
          store.path = plot.EQCC.dir,
          plot.CI = TRUE, save.separate.panels = FALSE)
# ----

## speed ratio ----
r.M2 <- list(EQCC = speedRat.M.d2$EQCC[qu.sel,],
             lower = speedRat.M.d2$lower[qu.sel,],
             upper = speedRat.M.d2$upper[qu.sel,])
r.M5 <- list(EQCC = speedRat.M.d5$EQCC[qu.sel,],
             lower = speedRat.M.d5$lower[qu.sel,],
             upper = speedRat.M.d5$upper[qu.sel,])
r.M11 <- list(EQCC = speedRat.M.d11$EQCC[qu.sel,],
              lower = speedRat.M.d11$lower[qu.sel,],
              upper = speedRat.M.d11$upper[qu.sel,])

r.F2 <- list(EQCC = speedRat.F.d2$EQCC[qu.sel,],
             lower = speedRat.F.d2$lower[qu.sel,],
             upper = speedRat.F.d2$upper[qu.sel,])
r.F5 <- list(EQCC = speedRat.F.d5$EQCC[qu.sel,],
             lower = speedRat.F.d5$lower[qu.sel,],
             upper = speedRat.F.d5$upper[qu.sel,])
r.F11 <- list(EQCC = speedRat.F.d11$EQCC[qu.sel,],
              lower = speedRat.F.d11$lower[qu.sel,],
              upper = speedRat.F.d11$upper[qu.sel,])

EQCC.plot(EQCC.out.list = list(r.F2,r.F5,r.F11,r.M2,r.M5,r.M11),
          y.var = "speedRat", qu.var = "time",
          y.label = "Speed variation",
          title.vec = c("EQCC speed variation (F: 10K)", "EQCC speed variation (F: Half)",
                        "EQCC speed variation (F: Full)", "EQCC speed variation (M: 10K)",
                        "EQCC speed variation (M: Half)", "EQCC speed variation (M: Full)",
                        "EQCC speed variation"),
          storeName.vec = c("speedRat-F-d2", "speedRat-F-d5", "speedRat-F-d11",
                            "speedRat-M-d2", "speedRat-M-d5", "speedRat-M-d11",
                            "speedRat-all"),
          box.labels = c("F: 10K", "F: Half", "F: Full",
                         "M: 10K", "M: Half", "M: Full"),
          box.locations = replicate(6,c(40,1.015), simplify = FALSE),
          qu.lvls = all.qu.lvls[qu.sel],
          dist.indices = 2:11,
          store.path = plot.EQCC.dir,
          plot.CI = TRUE, save.separate.panels = FALSE)
# ----

## EMQC curves: similar approach, using y.var = "ranking" and minor modifications
# to the plotting function

# ------------------------------------------------------------------------------

################################################################################



#### Part 2: Queensboro (QB) analysis ##########################################


## Data preparation ------------------------------------------------------------

# We work with a slightly modified population for the QB analysis, using some
# extra functions to redefine the data.


## Restriction to people with all interval speeds above certain threshold only
# currently used threshold: 8 km/h (valid VO2 formula => non-walking, too)
SomeWalking <- function(data){
  if (min(data$intSpeed) >= 8){
    return(rep(0, nrow(data)))
  } else { # walking on at least one of the other segments
    return(rep(1, nrow(data)))
  }
}

## alternative to the (more recent) allQuantiles() function
# requires proper definition of distance_vec
allQuantiles.old <- function(data, var_name){
  data$qu.lvl <- rep(0, dim(data)[1])
  data$qu.short <- rep(0, dim(data)[1])
  for (dist_index in 1:length(distance_vec)){
    current_dist_data <- data[near(data$dis, distance_vec[dist_index]), ]
    current_dist_ecdf <- ecdf(x = current_dist_data[[var_name]])
    current_qu_long <- current_dist_ecdf(current_dist_data[[var_name]])
    data[data$dis == distance_vec[dist_index], ]$qu.lvl <- current_qu_long
    
    # also append rounded version, that can be used for easier stratification
    # ensure that 0.020...01 up to 0.03 get assigned to 0.03 -> avoid 'round'
    current_qu_short <- 0.01*ceiling(100*current_qu_long)
    data[data$dis == distance_vec[dist_index], ]$qu.short <- current_qu_short  
  }
  
  return(data)
}

## QB quantifying metric: QBIntRat = Queensborough interval speed ratio v{4-6}/v{4}
QueensIntRatio <- function(data){
  data$intSpeed[5]/data$intSpeed[4]
}


## include only 9 distance points, omitting the irregular d5 and d8
distance_vec <- c(5,10,15,20,25,30,35,40,42.195)
NYC.long = read_csv('NYClong.csv')
NYC.long <- NYC.long[!NYC.long$distance %in% c("Half", "M20"),]
NYC.long <- mutate(NYC.long,
                   Sex = factor(Sex),
                   AgeGroup = factor(AgeGroup,
                                     levels = c("[18,19]","[20,29]","[30,39]",
                                                "[40,49]","[50,59]","[60,69]",
                                                "[70,79]", "[80,89]")),
                   distance = factor(distance))
NYC.long$Age2 = (NYC.long$Age - 18)/10
NYC.long$Mtime = rep(NYC.long$time[NYC.long$distance=="MAR"],each=9)

NYC.long$intSpeed = unlist(by(NYC.long,NYC.long$bib,IntervalSpeed))
NYC.long$cumSpeed = unlist(by(NYC.long,NYC.long$bib,CumulativeSpeed))
NYC.long$splitVec = unlist(by(NYC.long,NYC.long$bib,SplittingVector))
NYC.long$speedRat = unlist(by(NYC.long,NYC.long$bib,SpeedRatio))


## exclude the walkers/where formula VO2 not so accurate
NYC.long$walking = unlist(by(NYC.long,NYC.long$bib,SomeWalking))

NYC.runM <- NYC.long[NYC.long$Sex == "M" & NYC.long$walking == 0,]
NYC.runF <- NYC.long[NYC.long$Sex == "F" & NYC.long$walking == 0,]


## append quantiles and QBIntRat metric
NYC.runM <- allQuantiles.old(NYC.runM, 'time')
NYC.runF <- allQuantiles.old(NYC.runF, 'time')

NYC.runM$QBIntRat = rep(unlist(by(NYC.runM, NYC.runM$bib, QueensIntRatio)),
                        each = 9)
NYC.runF$QBIntRat = rep(unlist(by(NYC.runF, NYC.runF$bib, QueensIntRatio)),
                        each = 9)

NYC.runM <- select(NYC.runM,
                   c("bib", "Sex", "AgeGroup", "Age2", "dis", "time",
                     "Mtime", "intSpeed", "cumSpeed", "splitVec", "speedRat",
                     "qu.lvl", "qu.short", "QBIntRat"))

NYC.runF <- select(NYC.runF,
                   c("bib", "Sex", "AgeGroup", "Age2", "dis", "time",
                     "Mtime", "intSpeed", "cumSpeed", "splitVec", "speedRat",
                     "qu.lvl", "qu.short", "QBIntRat"))
# ------------------------------------------------------------------------------


## Trajectory visualisation ----------------------------------------------------
NYC.info <- data.frame(distance = c(20,21,22,23.5,25),
                       altitude = c(5, 19, 4, 9, 44))

ggplot(NYC.info) +
  geom_line(aes(x = distance, y = altitude), lwd = 1.5) +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)) +
  labs(x = "distance (km)",
       y = "altitude (m)") +
  geom_ribbon(aes(x = distance, ymin = 0, ymax = altitude), fill = "grey70")

ggsave(filename = paste0("QB-elevation.png"), path = plot.QB.dir,
       width = 6, height = 4, device='png')
# ------------------------------------------------------------------------------


## Proportion of runners with more slowing than expected (& crit. value) -------
# Ratio QBIntRat that one would expect physiologically, if one considers
# v_{15-20} as a flat segment, s.t. v_{4,f} = v_4 observed
QB.crit.flat <- (1/(1 +4.5*(14/sqrt(1000^2-14^2)*0.2 -0.2*15/
                              sqrt(1000^2-15^2) +
                              0.3*5/sqrt(1500^2-5^2) +
                              0.3*35/sqrt(1500^2-35^2))))

# proportion of runners with more QB slowing than expected
male.ecdf <- ecdf(NYC.runM$QBIntRat[near(NYC.runM$dis,42.195)])
female.ecdf <- ecdf(NYC.runF$QBIntRat[near(NYC.runF$dis,42.195)])

# male.ecdf(QB.crit.flat) # 0.6982
# female.ecdf(QB.crit.flat) # 0.5776


# aggregating 3 percentiles, and using v_4 as proxy for flat segment
# (k = 1  ~ 0.01, 0.02 and 0.03; k = 2 ~ 0.04, 0.05 and 0.06, etc.)
# note: nM and nF columns are just for verification purposes, not actually used
df.props.3.flat <- data.frame(
  qu = 0.03*(1:33),
  nM = sapply(1:33, function(k){
    nrow(NYC.runM[(0.03*k - 0.015 < NYC.runM$qu.short &
                     NYC.runM$qu.short < 0.03*k + 0.015) &
                    (near(NYC.runM$dis, 42.195)),])}),
  nF = sapply(1:33, function(k){
    nrow(NYC.runF[(0.03*k - 0.015 < NYC.runF$qu.short &
                     NYC.runF$qu.short < 0.03*k + 0.015) &
                    (near(NYC.runF$dis, 42.195)),])}),
  propM = sapply(1:33, function(k){
    df.aux <- NYC.runM[(0.03*k - 0.015 < NYC.runM$qu.short &
                          NYC.runM$qu.short < 0.03*k + 0.015) &
                         (near(NYC.runM$dis, 42.195)),]
    ecdf.aux <- ecdf(df.aux$QBIntRat)
    return(ecdf.aux(QB.crit.flat))
  }),
  propF = sapply(1:33, function(k){
    df.aux <- NYC.runF[(0.03*k - 0.015 < NYC.runF$qu.short &
                          NYC.runF$qu.short < 0.03*k + 0.015) &
                         (near(NYC.runF$dis, 42.195)),]
    ecdf.aux <- ecdf(df.aux$QBIntRat)
    return(ecdf.aux(QB.crit.flat))
  }))

## Bootstrapping to also build confidence intervals for these proportions
df.empty <- data.frame(matrix(NA, nrow = 33, ncol = 5))
colnames(df.empty) <- c("qu", "nM", "nF", "propM", "propF")
df.list <- replicate(100, df.empty, simplify = FALSE)


ids.M <- unique(NYC.runM$bib)
ids.F <- unique(NYC.runF$bib)

boot.rep <- 100
for (boot.index in 1:boot.rep){
  print(boot.index)
  set.seed(boot.index+123)
  new.ids.M <- sample(ids.M, size = length(ids.M), replace = TRUE)
  new.ids.F <- sample(ids.F, size = length(ids.F), replace = TRUE)
  
  # avoid logical indexing, where items selected multiple times appear only once
  NYC.resamM <- inner_join(x = NYC.runM, y = data.frame(bib = new.ids.M),
                           relationship = "many-to-many", by = join_by(bib))
  NYC.resamF <- inner_join(x = NYC.runF, y = data.frame(bib = new.ids.F),
                           relationship = "many-to-many", by = join_by(bib))
  
  # update quantile columns
  NYC.resamM <- allQuantiles.old(NYC.resamM, 'time')
  NYC.resamF <- allQuantiles.old(NYC.resamF, 'time')
  
  df.list[[boot.index]] <- data.frame(
    qu = 0.03*(1:33),
    nM = sapply(1:33, function(k){
      nrow(NYC.resamM[(0.03*k - 0.015 < NYC.resamM$qu.short &
                         NYC.resamM$qu.short < 0.03*k + 0.015) &
                        (near(NYC.resamM$dis, 42.195)),])}),
    nF = sapply(1:33, function(k){
      nrow(NYC.resamF[(0.03*k - 0.015 < NYC.resamF$qu.short &
                         NYC.resamF$qu.short < 0.03*k + 0.015) &
                        (near(NYC.resamF$dis, 42.195)),])}),
    propM = sapply(1:33, function(k){
      df.aux <- NYC.resamM[(0.03*k - 0.015 < NYC.resamM$qu.short &
                              NYC.resamM$qu.short < 0.03*k + 0.015) &
                             (near(NYC.resamM$dis, 42.195)),]
      ecdf.aux <- ecdf(df.aux$QBIntRat)
      return(ecdf.aux(QB.crit.flat))
    }),
    propF = sapply(1:33, function(k){
      df.aux <- NYC.resamF[(0.03*k - 0.015 < NYC.resamF$qu.short &
                              NYC.resamF$qu.short < 0.03*k + 0.015) &
                             (near(NYC.resamF$dis, 42.195)),]
      ecdf.aux <- ecdf(df.aux$QBIntRat)
      return(ecdf.aux(QB.crit.flat))
    }))
}


# Compute lower and upper bounds corresponding to the given significance;
# determine whether or not Bonferroni for multiple testing is to be applied

sgnf <- 0.05
propQB.title <- "propQB-flat-CI-3-aggreg.png"

# initialisation
df.propBounds <- data.frame(matrix(NA, nrow = 33, ncol = 4))
colnames(df.propBounds) <- c("lowerM", "lowerF", "upperM", "upperF")

for (qu.index in 1:33){
  # combine all (boot.rep) computed values at this entry, into vectors
  valsM <- sapply(seq_len(boot.rep), function(boot.index){
    return(df.list[[boot.index]]$propM[qu.index])
  })
  valsF <- sapply(seq_len(boot.rep), function(boot.index){
    return(df.list[[boot.index]]$propF[qu.index])
  })
  valsM.sorted <- sort(valsM)
  valsF.sorted <- sort(valsF)
  
  propM.lower <- valsM.sorted[max(floor(boot.rep*sgnf/2),1)]
  propF.lower <- valsF.sorted[max(floor(boot.rep*sgnf/2),1)]
  
  propM.upper <- valsM.sorted[max(floor(boot.rep*(1-sgnf/2)),1)]
  propF.upper <- valsF.sorted[max(floor(boot.rep*(1-sgnf/2)),1)]
  
  df.propBounds[qu.index,] <- c(propM.lower, propF.lower, propM.upper, propF.upper)
}

df.rearranged.3.flat <- data.frame(
  qu = rep(df.props.3.flat$qu,2),
  n = c(df.props.3.flat$nM, df.props.3.flat$nF),
  prop = c(df.props.3.flat$propM, df.props.3.flat$propF),
  propLower = c(df.propBounds$lowerM, df.propBounds$lowerF),
  propUpper = c(df.propBounds$upperM, df.propBounds$upperF),
  Sex = rep(c("M","F"), each = 33))


ggplot(df.rearranged.3.flat) +
  geom_line(aes(x = qu, y = prop, colour = Sex),
            lwd = 1.5) +
  geom_ribbon(aes(x = qu, ymin = propLower, ymax = propUpper,
                  fill = Sex, color = Sex),
              alpha = 0.3) +
  scale_colour_manual(name = "Sex", values = c("#41B6C4", "#0C2CB4")) +
  scale_fill_manual(name = "Sex", values = c("#41B6C4", "#0C2CB4")) +
  labs(x = "Finish time quantile",
       y = "") +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)) +
  plot_annotation(title = TeX(r'(Proportion with QBIntRat $< \kappa_{QB}$)'),
                  theme = theme(plot.title = element_text(hjust = 0.47),
                                text = element_text(size=15)))

ggsave(filename = propQB.title, path = plot.QB.dir,
       width = 6, height = 4, device='png')
# ------------------------------------------------------------------------------


## Introduction of Group 1 and Group 2, and QBIntRat exploration ---------------

# for separate analysis in Group 1 and Group 2 below.
NYC.M1 <- NYC.runM[NYC.runM$Mtime <= 195, ]
NYC.M2 <- NYC.runM[(NYC.runM$Mtime > 195 & NYC.runM$Mtime <= 265), ]

NYC.F1 <- NYC.runF[NYC.runF$Mtime <= 220, ]
NYC.F2 <- NYC.runF[(NYC.runF$Mtime > 220 & NYC.runF$Mtime <= 280), ]

# omit qu.lvl and qu.short columns, for application of EQCC functions
NYC.F1 <- select(NYC.F1, -c("qu.lvl", "qu.short"))
NYC.F2 <- select(NYC.F2, -c("qu.lvl", "qu.short"))

NYC.M1 <- select(NYC.M1, -c("qu.lvl", "qu.short"))
NYC.M2 <- select(NYC.M2, -c("qu.lvl", "qu.short"))


M1.ecdf <- ecdf(NYC.M1$QBIntRat[near(NYC.M1$dis,42.195)])
# M1.ecdf(QB.crit.flat) # 0.6080
M2.ecdf <- ecdf(NYC.M2$QBIntRat[near(NYC.M2$dis,42.195)])
# M2.ecdf(QB.crit.flat) # 0.7241

F1.ecdf <- ecdf(NYC.F1$QBIntRat[near(NYC.F1$dis,42.195)])
# F1.ecdf(QB.crit.flat) # 0.5387
F2.ecdf <- ecdf(NYC.F2$QBIntRat[near(NYC.F2$dis,42.195)])
# F2.ecdf(QB.crit.flat) # 0.5967

qu.sel.flat <- c(0.03,0.1,0.25,0.4,0.5,0.6,0.75,0.9,0.97)


# quantile values of QBIntRat - in qu.sel.flat
QBIntRat.quM1.flat <- round(quantile(NYC.M1$QBIntRat,
                                     probs = qu.sel.flat), digits = 4)
QBIntRat.quM2.flat <- round(quantile(NYC.M2$QBIntRat,
                                     probs = qu.sel.flat), digits = 4)

QBIntRat.quF1.flat <- round(quantile(NYC.F1$QBIntRat,
                                     probs = qu.sel.flat), digits = 4)
QBIntRat.quF2.flat <- round(quantile(NYC.F2$QBIntRat,
                                     probs = qu.sel.flat), digits = 4)

# xtable(rbind(QBIntRat.quM1.flat,
#              QBIntRat.quF1.flat,
#              QBIntRat.quM2.flat,
#              QBIntRat.quF2.flat), digits = 4)
# ------------------------------------------------------------------------------


## EQCC curves (with CI), conditioning on QBIntRat -----------------------------
# NOTE: aggregation by 3 percentiles! Cf. infra.

# slightly modified version of get.EQCC: aggregates 3 percentiles (to increase
# robustness in the case of small sample size) + modified distance defaults
get.EQCC.3 <- function(df, y.var, qu.var, dist.var = "dis",
                       qu.lvls = c(0.25,0.5,0.75), dist.indices = (1:9),
                       dist.vec = c(5,10,15,20,25,30,35,40,42.195),
                       cond.dist = 9){
  
  df.ext <- allQuantiles(data = df[
    near_in_vec(df[[dist.var]], dist.vec[dist.indices]),],
    var.name = qu.var, dist.var = dist.var, dist.vec = dist.vec,
    dist.indices = dist.indices, det.dist = cond.dist)
  
  df.EQCC <- data.frame(matrix(NA, nrow = length(qu.lvls),
                               ncol = length(dist.indices)))
  rownames(df.EQCC) <- paste0("qu.", qu.lvls)
  colnames(df.EQCC) <- paste0("D", dist.indices)
  
  for (qu.index in seq_along(qu.lvls)){
    # this.qu - 0.01, this.qu and this.qu + 0.01 will be included
    this.qu = df.ext[((df.ext$qu.det > (qu.lvls[qu.index] - 0.015)) &
                        (df.ext$qu.det < (qu.lvls[qu.index] + 0.015))),]
    
    if (nrow(this.qu) > 0){
      for (index in seq_along(dist.indices)){
        df.EQCC[qu.index, index] <- mean(this.qu[[y.var]][
          near(this.qu[[dist.var]], dist.vec[dist.indices[index]])])
      }
    }
  }
  
  return(df.EQCC)
}

# makes use of this modified get.EQCC.3 function, and selection of intSpeed and
# cumSpeed only (vs. more general multi-function) & other defaults for the
# distances and quantiles
EQCC.multi.3 <- function(df, qu.var, dist.var = "dis", id.var = "bib",
                         qu.lvls = c(0.03,0.05,0.1,0.2,0.25,0.3,0.4,0.5,
                                     0.6,0.7,0.75,0.8,0.9,0.95,0.97),
                         dist.vec = c(5,10,15,20,25,30,35,40,42.195),
                         cond.dist = 9,
                         comp.CI = FALSE, store.path = NA, store.name = NA,
                         sgnf = 0.05, boot.rep = 1000){
  
  ### Computation of EQCC values
  df.EQCC.intSpeed <- get.EQCC.3(df = df, y.var = "intSpeed", qu.var = qu.var,
                                 dist.var = dist.var, qu.lvls = qu.lvls,
                                 dist.indices = 1:9, dist.vec = dist.vec,
                                 cond.dist = cond.dist)
  df.EQCC.cumSpeed <- get.EQCC.3(df = df, y.var = "cumSpeed", qu.var = qu.var,
                                 dist.var = dist.var, qu.lvls = qu.lvls,
                                 dist.indices = 1:9, dist.vec = dist.vec,
                                 cond.dist = cond.dist)
  
  
  ### Computation of EQCC confidence intervals
  if (comp.CI){
    
    ### Compute EQCC values for boot.rep bootstrap resamples of the data
    
    df.empty <- data.frame(matrix(NA, nrow = length(qu.lvls), ncol = 9))
    rownames(df.empty) <- paste0("qu.", qu.lvls)
    colnames(df.empty) <- paste0("D", 1:9)
    
    df.list.intSpeed <- replicate(boot.rep, df.empty, simplify = FALSE)
    df.list.cumSpeed <- replicate(boot.rep, df.empty, simplify = FALSE)
    
    all.ids <- unique(df[[id.var]])
    
    for (boot.index in seq_len(boot.rep)){
      
      ## printing progress
      if (boot.index %% 25 == 0){
        print(boot.index)
      }
      
      ## bootstrap resample
      set.seed(boot.index)
      resample.ids <- sample(all.ids, size = length(all.ids), replace = TRUE)
      
      df.bibs <- data.frame(matrix(resample.ids, ncol = 1))
      colnames(df.bibs) <- id.var
      
      if (id.var == "bib"){ # avoid 'join_by' message by explicitly including it
        df.resampled <- inner_join(x = df, y = df.bibs,
                                   relationship = "many-to-many", by = join_by(bib))
      } else { # also works for id.var = "bib", but throws a message each time
        df.resampled <- inner_join(x = df, y = df.bibs,
                                   relationship = "many-to-many")
      }
      
      df.list.intSpeed[[boot.index]] <- get.EQCC.3(
        df = df.resampled, y.var = "intSpeed", qu.var = qu.var, dist.var = dist.var,
        qu.lvls = qu.lvls, dist.indices = 1:9, dist.vec = dist.vec,
        cond.dist = cond.dist)
      df.list.cumSpeed[[boot.index]] <- get.EQCC.3(
        df = df.resampled, y.var = "cumSpeed", qu.var = qu.var, dist.var = dist.var,
        qu.lvls = qu.lvls, dist.indices = 1:9, dist.vec = dist.vec,
        cond.dist = cond.dist)
    }
    
    ## Compute lower and upper bounds corresponding to the given significance
    df.lower.intSpeed <- df.empty
    df.upper.intSpeed <- df.empty

    df.lower.cumSpeed <- df.empty
    df.upper.cumSpeed <- df.empty
    
    for (qu.index in seq_along(qu.lvls)){
      for (dist.index in (1:9)){
        
        vals.intSpeed <- sapply(seq_len(boot.rep), function(boot.index){
          return(df.list.intSpeed[[boot.index]][qu.index,dist.index])
        })
        vals.cumSpeed <- sapply(seq_len(boot.rep), function(boot.index){
          return(df.list.cumSpeed[[boot.index]][qu.index,dist.index])
        })
        
        
        if (any(is.na(vals.intSpeed))){
          
          vals.intSpeed.noNA <- vals.intSpeed[!is.na(vals.intSpeed)]
          vals.intSpeed.sorted <- sort(vals.intSpeed.noNA)
          
          vals.cumSpeed.noNA <- vals.cumSpeed[!is.na(vals.cumSpeed)]
          vals.cumSpeed.sorted <- sort(vals.cumSpeed.noNA)
          
          rep.remaining <- length(vals.intSpeed.noNA)
          
          if (dist.index == 1){
            if (rep.remaining == (boot.rep - 1)){
              warning(paste0("At level qu = ", qu.lvls[qu.index],
                             ", 1 out of ", boot.rep, " bootstrap resamples",
                             " was omitted due to NA values."))
            } else {
              warning(paste0("At level qu = ", qu.lvls[qu.index], ", ",
                             boot.rep - rep.remaining,
                             " out of ", boot.rep, " bootstrap resamples",
                             " were omitted due to NA values."))
            }
          }
        } else {
          vals.intSpeed.sorted <- sort(vals.intSpeed)
          vals.cumSpeed.sorted <- sort(vals.cumSpeed)
          rep.remaining <- boot.rep
        }
        
        
        df.lower.intSpeed[qu.index, dist.index] <-
          vals.intSpeed.sorted[max(floor(rep.remaining*sgnf/2),1)]
        df.upper.intSpeed[qu.index, dist.index] <-
          vals.intSpeed.sorted[ceiling(rep.remaining*(1-sgnf/2))]
        
        df.lower.cumSpeed[qu.index, dist.index] <-
          vals.cumSpeed.sorted[max(floor(rep.remaining*sgnf/2),1)]
        df.upper.cumSpeed[qu.index, dist.index] <-
          vals.cumSpeed.sorted[ceiling(rep.remaining*(1-sgnf/2))]
      }
    }
    
  } else {
    df.lower.intSpeed <- NA
    df.upper.intSpeed <- NA
    
    df.lower.cumSpeed <- NA
    df.upper.cumSpeed <- NA
  }
  
  
  full.output <- list(
    resIntSpeed = list(EQCC = df.EQCC.intSpeed,
                       lower = df.lower.intSpeed, upper = df.upper.intSpeed),
    resCumSpeed = list(EQCC = df.EQCC.cumSpeed,
                       lower = df.lower.cumSpeed, upper = df.upper.cumSpeed)
  )
  
  ### Storage of results
  if (!is.na(store.path)){
    
    if (! dir.exists(store.path)){
      dir.create(store.path)
    }
    
    list.save(full.output,
              file = paste0(store.path, "/full.output.", store.name, ".Rdata"))
  }
  
  
  ### Also return results to console
  return(full.output)
}


all.qu.lvls.flat <- c(0.03,0.05,0.1,0.2,0.25,0.3,0.4,0.45,0.5,0.55,
                      0.6,0.7,0.75,0.8,0.9,0.95,0.97)
# don't look below 0.02 or above 0.99: also one below/above should be present
# with aggregation of 3 percentiles
qu.sel <- c(1,3,5,7,9,11,13,15,17)



# cond.dist = random index, as QBIntRat is distance-independent.
EQCC.M1.intRat <- EQCC.multi.3(df = NYC.M1, qu.var = "QBIntRat",
                               qu.lvls = all.qu.lvls.flat,
                               cond.dist = 9, comp.CI = TRUE,
                               store.path = store.QB.EQCC.dir,
                               store.name = "M-gr1",
                               boot.rep = 1000)

EQCC.M2.intRat <- EQCC.multi.3(df = NYC.M2, qu.var = "QBIntRat",
                               qu.lvls = all.qu.lvls.flat,
                               cond.dist = 9, comp.CI = TRUE,
                               store.path = store.QB.EQCC.dir,
                               store.name = "M-gr2",
                               boot.rep = 1000)

EQCC.F1.intRat <- EQCC.multi.3(df = NYC.F1, qu.var = "QBIntRat",
                               qu.lvls = all.qu.lvls.flat,
                               cond.dist = 9, comp.CI = TRUE,
                               store.path = store.QB.EQCC.dir,
                               store.name = "F-gr1",
                               boot.rep = 1000)

EQCC.F2.intRat <- EQCC.multi.3(df = NYC.F2, qu.var = "QBIntRat",
                               qu.lvls = all.qu.lvls.flat,
                               cond.dist = 9, comp.CI = TRUE,
                               store.path = store.QB.EQCC.dir,
                               store.name = "F-gr2",
                               boot.rep = 1000)

intSpeed.M1.intRat <- EQCC.M1.intRat$resIntSpeed
intSpeed.F1.intRat <- EQCC.F1.intRat$resIntSpeed
intSpeed.M2.intRat <- EQCC.M2.intRat$resIntSpeed
intSpeed.F2.intRat <- EQCC.F2.intRat$resIntSpeed

cumSpeed.M1.intRat <- EQCC.M1.intRat$resCumSpeed
cumSpeed.F1.intRat <- EQCC.F1.intRat$resCumSpeed
cumSpeed.M2.intRat <- EQCC.M2.intRat$resCumSpeed
cumSpeed.F2.intRat <- EQCC.F2.intRat$resCumSpeed




## interval speed
v.intRat.M1 <- list(EQCC = intSpeed.M1.intRat$EQCC[qu.sel,],
                    lower = intSpeed.M1.intRat$lower[qu.sel,],
                    upper = intSpeed.M1.intRat$upper[qu.sel,])
v.intRat.F1 <- list(EQCC = intSpeed.F1.intRat$EQCC[qu.sel,],
                    lower = intSpeed.F1.intRat$lower[qu.sel,],
                    upper = intSpeed.F1.intRat$upper[qu.sel,])

v.intRat.M2 <- list(EQCC = intSpeed.M2.intRat$EQCC[qu.sel,],
                    lower = intSpeed.M2.intRat$lower[qu.sel,],
                    upper = intSpeed.M2.intRat$upper[qu.sel,])

v.intRat.F2 <- list(EQCC = intSpeed.F2.intRat$EQCC[qu.sel,],
                    lower = intSpeed.F2.intRat$lower[qu.sel,],
                    upper = intSpeed.F2.intRat$upper[qu.sel,])


# Separate plot for group 1 and group 2
EQCC.plot(EQCC.out.list = list(v.intRat.F1,v.intRat.M1),
          y.var = "intSpeed", qu.var = "QBIntRat",
          y.label = "Interval speed (km/h)",
          title.vec = c("EQCC interval speed (F: Group 1)",
                        "EQCC interval speed (M: Group 1)",
                        "EQCC interval speed conditional on QBIntRat"),
          storeName.vec = c("intSpeed-F-gr1",
                            "intSpeed-M-gr1",
                            "intSpeed-all-gr1"),
          box.labels = c("F: Group 1","M: Group 1"),
          box.locations = replicate(6,c(38,16.1), simplify = FALSE), # old/mix
          # box.locations = replicate(6,c(40,13.7), simplify = FALSE), # quartiles
          qu.lvls = all.qu.lvls.flat[qu.sel],
          dist.indices = 1:9,
          dist.vec = c(5,10,15,20,25,30,35,40,42.195),
          store.path = plot.QB.EQCC.dir,
          plot.CI = TRUE, save.separate.panels = FALSE)
EQCC.plot(EQCC.out.list = list(v.intRat.F2,v.intRat.M2),
          y.var = "intSpeed", qu.var = "QBIntRat",
          y.label = "Interval speed (km/h)",
          title.vec = c("EQCC interval speed (F: Group 2)",
                        "EQCC interval speed (M: Group 2)",
                        "EQCC interval speed conditional on QBIntRat"),
          storeName.vec = c("intSpeed-F-gr2",
                            "intSpeed-M-gr2",
                            "intSpeed-all-gr2"),
          box.labels = c("F: Group 2","M: Group 2"),
          box.locations = replicate(6,c(38,13), simplify = FALSE), # old/mix
          # box.locations = replicate(6,c(40,13.7), simplify = FALSE), # quartiles
          qu.lvls = all.qu.lvls.flat[qu.sel],
          dist.indices = 1:9,
          dist.vec = c(5,10,15,20,25,30,35,40,42.195),
          store.path = plot.QB.EQCC.dir,
          plot.CI = TRUE, save.separate.panels = FALSE)




## cumulative speed
c.intRat.M1 <- list(EQCC = cumSpeed.M1.intRat$EQCC[qu.sel,],
                    lower = cumSpeed.M1.intRat$lower[qu.sel,],
                    upper = cumSpeed.M1.intRat$upper[qu.sel,])
c.intRat.F1 <- list(EQCC = cumSpeed.F1.intRat$EQCC[qu.sel,],
                    lower = cumSpeed.F1.intRat$lower[qu.sel,],
                    upper = cumSpeed.F1.intRat$upper[qu.sel,])

c.intRat.M2 <- list(EQCC = cumSpeed.M2.intRat$EQCC[qu.sel,],
                    lower = cumSpeed.M2.intRat$lower[qu.sel,],
                    upper = cumSpeed.M2.intRat$upper[qu.sel,])
c.intRat.F2 <- list(EQCC = cumSpeed.F2.intRat$EQCC[qu.sel,],
                    lower = cumSpeed.F2.intRat$lower[qu.sel,],
                    upper = cumSpeed.F2.intRat$upper[qu.sel,])

# Plotting could be performed for cumSpeed, too, but we look only at the overall
# cumulative speed at full marathon distance here, that can be matched to the
# intSpeed EQCC curves.


# xtable(rbind(c.intRat.M1$lower[,9],
#              c.intRat.M1$EQCC[,9],
#              c.intRat.M1$upper[,9],
#              c.intRat.F1$lower[,9],
#              c.intRat.F1$EQCC[,9],
#              c.intRat.F1$upper[,9]), digits = 2)

# xtable(rbind(c.intRat.M2$lower[,9],
#              c.intRat.M2$EQCC[,9],
#              c.intRat.M2$upper[,9],
#              c.intRat.F2$lower[,9],
#              c.intRat.F2$EQCC[,9],
#              c.intRat.F2$upper[,9]), digits = 2)

# ------------------------------------------------------------------------------

################################################################################



#### Part 3: Sensitivity analysis to the QB analysis ###########################

# For the sensitivity analysis reported in the paper, rerun the above QB analysis
# with an alternative definition of the population, i.e. other walking threshold,
# and other subdivision in Group 1, Group 2 (and remainder)

# NOTE: make sure that the data to start with (NYC.long) are as in Part 2, i.e.,
# omitting D5 and D8, with distance_vec of length 9 only

## Data preparation ------------------------------------------------------------

SomeWalkingSens <- function(data){
  if (min(data$intSpeed) >= 7.5){
    return(rep(0, nrow(data)))
  } else {
    return(rep(1, nrow(data)))
  }
}

NYC.long$walkingSens = unlist(by(NYC.long,NYC.long$bib,SomeWalkingSens))

NYC.runM.sens <- NYC.long[NYC.long$Sex == "M" & NYC.long$walkingSens == 0,]
NYC.runF.sens <- NYC.long[NYC.long$Sex == "F" & NYC.long$walkingSens == 0,]


## append quantiles and QBIntRat
NYC.runM.sens <- allQuantiles.old(NYC.runM.sens, 'time')
NYC.runF.sens <- allQuantiles.old(NYC.runF.sens, 'time')

NYC.runM.sens$QBIntRat = rep(unlist(by(NYC.runM.sens, NYC.runM.sens$bib, QueensIntRatio)),
                            each = 9)
NYC.runF.sens$QBIntRat = rep(unlist(by(NYC.runF.sens, NYC.runF.sens$bib, QueensIntRatio)),
                            each = 9)

NYC.runM.sens <- select(NYC.runM.sens,
                       c("bib", "Sex", "AgeGroup", "Age2", "dis", "time",
                         "Mtime", "intSpeed", "cumSpeed", "splitVec", "speedRat",
                         "qu.lvl", "qu.short", "QBIntRat"))

NYC.runF.sens <- select(NYC.runF.sens,
                       c("bib", "Sex", "AgeGroup", "Age2", "dis", "time",
                         "Mtime", "intSpeed", "cumSpeed", "splitVec", "speedRat",
                         "qu.lvl", "qu.short", "QBIntRat"))


NYC.sensM1 <- NYC.runM.sens[NYC.runM.sens$Mtime <= 200, ]
NYC.sensM2 <- NYC.runM.sens[(NYC.runM.sens$Mtime > 200 & NYC.runM.sens$Mtime <= 280), ]

NYC.sensF1 <- NYC.runF.sens[NYC.runF.sens$Mtime <= 225, ]
NYC.sensF2 <- NYC.runF.sens[(NYC.runF.sens$Mtime > 225 & NYC.runF.sens$Mtime <= 290), ]

# make sure that proportions of male and female within group 1 and 2 stay roughly
# equal, otherwise we are comparing a different population for male vs female

# nrow(NYC.M1)/9
# nrow(NYC.sensM1)/9
# 
# nrow(NYC.F1)/9
# nrow(NYC.sensF1)/9
# 
# nrow(NYC.M2)/9
# nrow(NYC.sensM2)/9
# 
# nrow(NYC.F2)/9
# nrow(NYC.sensF2)/9

# nrow(NYC.M1)/nrow(NYC.M2)
# nrow(NYC.F1)/nrow(NYC.F2)
# 
# nrow(NYC.sensM1)/nrow(NYC.sensM2)
# nrow(NYC.sensF1)/nrow(NYC.sensF2)
# ------------------------------------------------------------------------------


## QBIntRat exploration --------------------------------------------------------
qu.sel.flat <- c(0.03,0.1,0.25,0.4,0.5,0.6,0.75,0.9,0.97)

# quantile values of QBIntRat - in qu.sel.flat
QBIntRat.quSensM1.flat <- round(quantile(NYC.sensM1$QBIntRat,
                                         probs = qu.sel.flat), digits = 4)
QBIntRat.quSensM2.flat <- round(quantile(NYC.sensM2$QBIntRat,
                                         probs = qu.sel.flat), digits = 4)

QBIntRat.quSensF1.flat <- round(quantile(NYC.sensF1$QBIntRat,
                                         probs = qu.sel.flat), digits = 4)
QBIntRat.quSensF2.flat <- round(quantile(NYC.sensF2$QBIntRat,
                                         probs = qu.sel.flat), digits = 4)

# xtable(rbind(QBIntRat.quSensM1.flat,
#              QBIntRat.quSensF1.flat,
#              QBIntRat.quSensM2.flat,
#              QBIntRat.quSensF2.flat), digits = 4)
# ------------------------------------------------------------------------------


## EQCC curves (with CI), conditioning on QBIntRat -----------------------------
EQCC.sensM1.intRat <- EQCC.multi.3(df = NYC.sensM1, qu.var = "QBIntRat",
                                  qu.lvls = all.qu.lvls.flat,
                                  cond.dist = 9, comp.CI = TRUE,
                                  store.path = store.QB.EQCC.sens.dir,
                                  store.name = "M-gr1-sens",
                                  boot.rep = 1000)

EQCC.sensM2.intRat <- EQCC.multi.3(df = NYC.sensM2, qu.var = "QBIntRat",
                                  qu.lvls = all.qu.lvls.flat,
                                  cond.dist = 9, comp.CI = TRUE,
                                  store.path = store.QB.EQCC.sens.dir,
                                  store.name = "M-gr2-sens",
                                  boot.rep = 1000)

EQCC.sensF1.intRat <- EQCC.multi.3(df = NYC.sensF1, qu.var = "QBIntRat",
                                  qu.lvls = all.qu.lvls.flat,
                                  cond.dist = 9, comp.CI = TRUE,
                                  store.path = store.QB.EQCC.sens.dir,
                                  store.name = "F-gr1-sens",
                                  boot.rep = 1000)

EQCC.sensF2.intRat <- EQCC.multi.3(df = NYC.sensF2, qu.var = "QBIntRat",
                                  qu.lvls = all.qu.lvls.flat,
                                  cond.dist = 9, comp.CI = TRUE,
                                  store.path = store.QB.EQCC.sens.dir,
                                  store.name = "F-gr2-sens",
                                  boot.rep = 1000)


intSpeed.sensM1.intRat <- EQCC.sensM1.intRat$resIntSpeed
intSpeed.sensF1.intRat <- EQCC.sensF1.intRat$resIntSpeed
intSpeed.sensM2.intRat <- EQCC.sensM2.intRat$resIntSpeed
intSpeed.sensF2.intRat <- EQCC.sensF2.intRat$resIntSpeed

cumSpeed.sensM1.intRat <- EQCC.sensM1.intRat$resCumSpeed
cumSpeed.sensF1.intRat <- EQCC.sensF1.intRat$resCumSpeed
cumSpeed.sensM2.intRat <- EQCC.sensM2.intRat$resCumSpeed
cumSpeed.sensF2.intRat <- EQCC.sensF2.intRat$resCumSpeed



qu.sel <- c(1,3,9,15,17)


## interval speed
v.intRat.sensM1 <- list(EQCC = intSpeed.sensM1.intRat$EQCC[qu.sel,],
                       lower = intSpeed.sensM1.intRat$lower[qu.sel,],
                       upper = intSpeed.sensM1.intRat$upper[qu.sel,])
v.intRat.sensF1 <- list(EQCC = intSpeed.sensF1.intRat$EQCC[qu.sel,],
                       lower = intSpeed.sensF1.intRat$lower[qu.sel,],
                       upper = intSpeed.sensF1.intRat$upper[qu.sel,])

v.intRat.sensM2 <- list(EQCC = intSpeed.sensM2.intRat$EQCC[qu.sel,],
                       lower = intSpeed.sensM2.intRat$lower[qu.sel,],
                       upper = intSpeed.sensM2.intRat$upper[qu.sel,])

v.intRat.sensF2 <- list(EQCC = intSpeed.sensF2.intRat$EQCC[qu.sel,],
                       lower = intSpeed.sensF2.intRat$lower[qu.sel,],
                       upper = intSpeed.sensF2.intRat$upper[qu.sel,])


# Separate plot for group 1 and group 2
EQCC.plot(EQCC.out.list = list(v.intRat.sensF1,v.intRat.sensM1),
          y.var = "intSpeed", qu.var = "QBIntRat",
          y.label = "Interval speed (km/h)",
          title.vec = c("EQCC interval speed (F: Group 1')",
                        "EQCC interval speed (M: Group 1')",
                        "EQCC interval speed conditional on QBIntRat"),
          storeName.vec = c("intSpeed-F-gr1-sens",
                            "intSpeed-M-gr1-sens",
                            "intSpeed-all-gr1-sens"),
          box.labels = c("F: Group 1'","M: Group 1'"),
          box.locations = replicate(6,c(38,15.9), simplify = FALSE),
          qu.lvls = all.qu.lvls.flat[qu.sel],
          dist.indices = 1:9,
          dist.vec = c(5,10,15,20,25,30,35,40,42.195),
          store.path = plot.QB.EQCC.sens.dir,
          plot.CI = TRUE, save.separate.panels = FALSE)
EQCC.plot(EQCC.out.list = list(v.intRat.sensF2,v.intRat.sensM2),
          y.var = "intSpeed", qu.var = "QBIntRat",
          y.label = "Interval speed (km/h)",
          title.vec = c("EQCC interval speed (F: Group 2')",
                        "EQCC interval speed (M: Group 2')",
                        "EQCC interval speed conditional on QBIntRat"),
          storeName.vec = c("intSpeed-F-gr2-sens",
                            "intSpeed-M-gr2-sens",
                            "intSpeed-all-gr2-sens"),
          box.labels = c("F: Group 2'","M: Group 2'"),
          box.locations = replicate(6,c(38,12.5), simplify = FALSE),
          qu.lvls = all.qu.lvls.flat[qu.sel],
          dist.indices = 1:9,
          dist.vec = c(5,10,15,20,25,30,35,40,42.195),
          store.path = plot.QB.EQCC.sens.dir,
          plot.CI = TRUE, save.separate.panels = FALSE)






## cumulative speed
c.intRat.sensM1 <- list(EQCC = cumSpeed.sensM1.intRat$EQCC[qu.sel,],
                       lower = cumSpeed.sensM1.intRat$lower[qu.sel,],
                       upper = cumSpeed.sensM1.intRat$upper[qu.sel,])
c.intRat.sensF1 <- list(EQCC = cumSpeed.sensF1.intRat$EQCC[qu.sel,],
                       lower = cumSpeed.sensF1.intRat$lower[qu.sel,],
                       upper = cumSpeed.sensF1.intRat$upper[qu.sel,])

c.intRat.sensM2 <- list(EQCC = cumSpeed.sensM2.intRat$EQCC[qu.sel,],
                       lower = cumSpeed.sensM2.intRat$lower[qu.sel,],
                       upper = cumSpeed.sensM2.intRat$upper[qu.sel,])
c.intRat.sensF2 <- list(EQCC = cumSpeed.sensF2.intRat$EQCC[qu.sel,],
                       lower = cumSpeed.sensF2.intRat$lower[qu.sel,],
                       upper = cumSpeed.sensF2.intRat$upper[qu.sel,])

# xtable(rbind(c.intRat.sensM1$lower[,9],
#              c.intRat.sensM1$EQCC[,9],
#              c.intRat.sensM1$upper[,9],
#              c.intRat.sensF1$lower[,9],
#              c.intRat.sensF1$EQCC[,9],
#              c.intRat.sensF1$upper[,9]))
# 
# xtable(rbind(c.intRat.sensM2$lower[,9],
#              c.intRat.sensM2$EQCC[,9],
#              c.intRat.sensM2$upper[,9],
#              c.intRat.sensF2$lower[,9],
#              c.intRat.sensF2$EQCC[,9],
#              c.intRat.sensF2$upper[,9]), digits = 2)
# ------------------------------------------------------------------------------

################################################################################




