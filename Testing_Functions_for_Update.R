# Testing out the functions:
devtools::load_all()

# phyloseq2table ---------------------------------------------------------------
phyloseq2table(physeq = physeq_S2D2,
               include_sample_data=TRUE,
               sample_col_keep = c('IS_CONTROL', 'Buoyant_density', "Microcosm_replicate", "Substrate"),
               control_var = "Substrate",
               control_expr = "12C-Con") # this works.

# qSIP_atom_excess_format ------------------------------------------------
head(sample_data(physeq_S2D2))

qSIP_atom_excess_format(physeq = physeq_S2D2,
                        control_var = "Substrate",
                        control_expr ="12C-Con",
                        treatment_rep = "Microcosm_replicate")




  # formatting input
# create new columns:
  cols <- c('IS_CONTROL', 'Buoyant_density', "Microcosm_replicate", "Substrate")
  control_expr ="12C-Con"
  # Create dataframe with Phyloseq2table and the expression for identifying the control
  df_OTU <- phyloseq2table(physeq_S2D2,
                          include_sample_data=TRUE,
                          sample_col_keep=cols,
                          control_var = "Substrate",
                          control_expr="12C-Con")

  # removing 'infinite' BD values
  tmp <- colnames(df_OTU)
  # Update the df_OTU table
  df_OTU <- df_OTU %>%
    dplyr::mutate(Buoyant_density = as.numeric(as.character(Buoyant_density)),
                  Count = as.numeric(as.character(Count))) %>%
    dplyr::filter(! is.infinite(Buoyant_density)) %>%
    dplyr::filter(! is.na(Buoyant_density)) %>%
    as.data.frame()
  colnames(df_OTU) = tmp

  # return
  return(df_OTU)

# qSIP_atom_excess -----------------------------------------------------------
  devtools::load_all()

  x <- qSIP_atom_excess(physeq = physeq_S2D2,
                        control_var = "Substrate",
                   control_expr =   "12C-Con",
                   treatment_rep = "Microcosm_replicate",
                   isotope='13C',
                   df_OTU_W=NULL)

  # For testing: Parameters
  physeq = physeq_S2D2
  control_expr =  "12C-Con"
  control_var = "Substrate"
  treatment_rep = "Microcosm_replicate"
  isotope='13C'
  df_OTU_W=NULL

    # formatting input
    if(is.null(df_OTU_W)){
      no_boot = TRUE
    } else {
      no_boot = FALSE
    }

    if(no_boot == TRUE){
      df_OTU <- qSIP_atom_excess_format(physeq,
                                        control_expr = control_expr,
                                        treatment_rep = treatment_rep,
                                        control_var = control_var)
      if(nrow(df_OTU) == 0){
        stop('No rows in OTU table after qSIP_atom_excess_format(). Check control_exp & treatment_rep')
      }

      # *Testing*:
      names(df_OTU)
      # colnames: OTU = Otu; SAMPLE_JOIN = sample names; IS_CONTROL = expression to determine controls

      # BD shift (Z) from the OTU table with the associated densities
      df_OTU_W <- df_OTU %>%
        # weighted mean buoyant density (W)
        dplyr::mutate(Buoyant_density = as.numeric(as.character(Buoyant_density)),
                       Count = as.numeric(as.character(Count))) %>%
        # IS_CONTROL should be a boolean vector? !
        dplyr::group_by(IS_CONTROL, OTU, .data[[treatment_rep]]) %>%
        # get mean buoyant density per control/Treatment and OTU over all replicates:
        dplyr::summarize(W = stats::weighted.mean(Buoyant_density, Count, na.rm=TRUE)) %>%
        dplyr::ungroup()
    }

    df_OTU_s <- df_OTU_W %>%
      # mean W of replicate gradients
      dplyr::group_by(IS_CONTROL, OTU) %>%
      dplyr::summarize(Wm = mean(W, na.rm=TRUE)) %>%
      # BD shift (Z)
      dplyr::group_by(OTU) %>%
      dplyr::mutate(IS_CONTROL = ifelse(IS_CONTROL==TRUE, 'Wlight', 'Wlab')) %>%
      tidyr::spread(IS_CONTROL, Wm) %>%
      dplyr::mutate(Z = Wlab - Wlight) %>%
      dplyr::ungroup()

    # atom excess (A)
    ## pt1: Calculate G+C content and Light Molecular Weight
    df_OTU_s <- df_OTU_s %>%
      dplyr::mutate(
        Gi = calc_Gi(Wlight),
        Mlight = 0.496 * Gi + 307.691
      )
    ## pt2: Calculate Max Heavy Molecular Weight
    df_OTU_s <- df_OTU_s %>%
      dplyr::mutate(
        Mheavymax = calc_Mheavymax(Mlight = Mlight, Gi = Gi, isotope = isotope)
      )
    ## pt3: Calculate Lab Molecular Weight and Atom Excess (A)
    df_OTU_s <- df_OTU_s %>%
      dplyr::mutate(
        # Calculate Mlab first so it's available for the next line
        Mlab = (Z / Wlight + 1) * Mlight,
        # Calculate Atom Excess (A)
        A = calc_atom_excess(Mlab = Mlab, Mlight = Mlight,
                             Mheavymax = Mheavymax, isotope = isotope)
      )

    ## flow control: bootstrap and return results
    if(no_boot){
      return(list(W=df_OTU_W, A=df_OTU_s))
    } else {
      return(df_OTU_s)


  # sampling with replacement from control & treatment for each OTU
  sample_W = function(df, n_sample){
    n_light = n_sample[1]
    n_lab = n_sample[2]
    # parsing df
    df_light = df[df$IS_CONTROL==TRUE,]
    df_lab = df[df$IS_CONTROL==FALSE,]
    # sampling
    if(length(df_light$W) > 1){
      W_light = base::sample(df_light$W, n_light, replace=TRUE)
    } else {
      W_light = rep(df_light$W, n_light)
    }
    if(length(df_lab$W) > 1){
      W_lab = base::sample(df_lab$W, n_lab, replace=TRUE)
    } else {
      W_lab = rep(df_lab$W, n_lab)
    }
    # creating new data.frames
    df_light = data.frame(IS_CONTROL=TRUE,
                          #OTU=rep(otu, n_light),
                          W = W_light)
    df_lab = data.frame(IS_CONTROL=FALSE,
                        #OTU=rep(otu, n_lab),
                        W = W_lab)
    return(rbind(df_light, df_lab))
  }


  # shuffling weighted mean densities (W)
  .qSIP_bootstrap = function(atomX,
                             isotope='13C',
                             n_sample=c(3,3),
                             bootstrap_id = 1){
    # making a new (subsampled with replacement) dataset
    n_sample = c(3,3)  # control, treatment
    dots = stats::setNames(list(~lapply(data, sample_W, n_sample=n_sample)), "ret")
    df_OTU_W = atomX$W %>%
      dplyr::group_by_("OTU") %>%
      tidyr::nest() %>%
      dplyr::mutate_(.dots=dots) %>%
      dplyr::select_("-data") %>%
      tidyr::unnest()

    # calculating atom excess
    atomX = qSIP_atom_excess(physeq=NULL,
                             df_OTU_W=df_OTU_W,
                             control_expr=NULL,
                             treatment_rep=NULL,
                             isotope=isotope)
    atomX$bootstrap_id = bootstrap_id
    return(atomX)
  }


# qSIP_bootstrap----------------------------------------------------------------
