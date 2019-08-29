library(shiny)
library(swCRTdesign)#Built under version 2.2/3.0
library(ggplot2)


#https://www.shinyapps.io
#rsconnect::deployApp('P:/Year_2/SWCRT RA/Shiny/Shiny_app_uploaded')
#rsconnect::deployApp('P:/Year_2/SWCRT RA/Shiny/stepped_wedge_power_calculation')
#hint: dog



#5/20 (_5): BIG SHIFT: decided to not depend on CRAN functions in the end (will have to add the helper functions here and remove library(swCRTdesign)!)
#This will allow me to mute annoying error messages
#Emily: move swDsn to package version; keep swPwr.FlexN
#Same code as v 3.0

####################################################
##Busy indicator 
####################################################

#Shinysky code, modified to be a little simpler and not display a gif
#Referenced: https://github.com/AnalytixWare/ShinySky
busyIndicator <- function(text = "Calculation in progress..",img = NA, wait=1000) {
  tagList(
    singleton(tags$head(
      tags$link(rel="stylesheet", type="text/css",href="shinysky/busyIndicator/busyIndicator.css")
    ))
    ,div(class="shinysky-busy-indicator",p(text),img(src=img))
    ,tags$script(sprintf(
      "	setInterval(function(){
  		 	 if ($('html').hasClass('shiny-busy')) {
  		    setTimeout(function() {
  		      if ($('html').hasClass('shiny-busy')) {
  		        $('div.shinysky-busy-indicator').show()
  		      }
  		    }, %d)  		    
  		  } else {
  		    $('div.shinysky-busy-indicator').hide()
  		  }
  		},100)
  		",wait)
    )
  )	
}

####################################################
##swDsn
####################################################

#Build the design matrix (see documentation for R package swCRTdesign)
#This should be identical to the package
swDsn <- function (clusters, tx.effect = 1, extra.time = 0, all.ctl.time0 = TRUE) 
{
  clusters.per.wave <- clusters[clusters != 0]#edit: clusters[clusters != 0], so this is really clusters per wave and not clusters per time point
  total.clusters <- sum(clusters.per.wave)
  waves <- length(clusters.per.wave)
  total.time <- (length(clusters) + 1 + extra.time)#edit: length(clusters) instead of waves
  swBlk <- round(upper.tri(matrix(1, length(clusters), total.time)))#edit: length(clusters) instead of waves
  swPreDesign.MAT <- cbind(clusters, swBlk)#edit: clusters instead of clusters.per.wave
  swPreDesign.LIST <- apply(swPreDesign.MAT, 1, function(z) {
    rep(z[-1], z[1])
  })
  swDsn <- matrix(unlist(swPreDesign.LIST), ncol = total.time, 
                  byrow = TRUE)
  if (!all.ctl.time0) {
    swBlk <- swBlk[, -1]
    swDsn <- swDsn[, -1]
    #clusters.per.wave <- clusters.per.wave[-c(waves)]
    #total.clusters <- sum(clusters.per.wave)
    #waves <- length(clusters.per.wave)
    total.time <- dim(swDsn)[2]
  }
  if (!is.null(tx.effect)) {
    swTxEffectPreDesign.LIST <- apply(swDsn, 1, function(z) {
      Nctl <- sum(z == 0)
      rowWithTxEffect <- c(rep(0, Nctl), tx.effect, rep(1, 
                                                        length(z)))[1:length(z)]
      rowWithTxEffect
    })
    swTxEffectDsn <- matrix(unlist(swTxEffectPreDesign.LIST), 
                            ncol = total.time, byrow = TRUE)
    swDsn <- swTxEffectDsn
    swBlk <- unique(swDsn)
  }
  list(swDsn = swDsn, swDsn.unique.clusters = swBlk, n.waves = waves, 
       clusters = clusters.per.wave, n.clusters = total.clusters, 
       tx.effect = tx.effect, total.time = total.time, extra.time = extra.time)
}

####################################################
##swPwr
####################################################

#Compared to the package version, this doesn't have as many warnings and the ICC/CAC isn't built in
#Body of the code is almost identical to swCRTdesign package v. 3.0; 
#biggest differences to look out for when making updates are in the inputs (different order and defaults) and warnings

swPwr.FlexN <- function (design, distn, n, mu0, mu1, tau, eta = 0, rho = 0, sigma, gamma = 0, 
                         alpha = 0.05, retDATA = FALSE) 
{
  obs.per.cluster.per.time <- n
  if (length(n)>1){
    warning("When sample sizes are not uniform, power depends on order of clusters (see documentation).")
  }
  if (rho != 0 & (tau == 0 | eta == 0)){
    stop("If either tau or eta is zero, the associated correlation (rho) must be zero.")
  }
  theta <- (mu1 - mu0)
  muBar <- (mu0 + mu1)/2
  if (distn == "gaussian"){ 
    sigSq <- sigma^2
  }else if (distn == "binomial") {
    if (!missing(sigma)) 
      warning("sigma is not used when distn=binomial")
    sigSq <- muBar * (1 - muBar)
    sigma <- NA
    if ((tau^2 + eta^2 + gamma^2) > sigSq) 
      stop("tau^2 + eta^2 + gamma^2 must be less than muBar*(1-muBar) when distn=binomial")
  }
  if (length(tau) > 1 | length(eta) > 1) 
    stop("Function cannot compute stepped wedge design Power for tau-vector or eta-vector; tau and eta must be scalars.")#EDIT: changed to stop from warning
  I.rep <- design$clusters
  I <- design$n.clusters
  J <- design$total.time
  swDesign <- design$swDsn
  swDesignUnique <- design$swDsn.unique.clusters
  if (any(rowSums(swDesign) == 0)) 
    warning("For the specified total number of clusters (I), total number of time periods (J), and number of cluster repetitions (I.rep), the specified stepped wedge design has at least one cluster which does not crossover from control(0) to treatment(1) arm.")
  X.ij <- as.vector(t(swDesign))
  beta.blk <- rbind(diag(1, J - 1, J - 1), 0)
  Xmat.blk <- matrix(rep(as.vector(t(cbind(1, beta.blk))), I), ncol = J, byrow = TRUE)
  Xmat <- cbind(Xmat.blk, X.ij)
  if (length(n) == 1){
    Wmat.blk <- tau^2 + diag(sigSq/n+gamma^2, J)#matrix V_i
    Wmat.partial <- kronecker(diag(1, I), Wmat.blk)#matrix V
    #above, I is the number of times that Wmat.blk gets repeated along the diagonal of a matrix filled with 0's otherwise
  }else{
    #Took this nMat creation from swSim: each row is for a specific cluster
    if ((is.vector(n) & length(n) > 1)) {
      #n is a vector with one entry per cluster
      if (length(n) != design$n.clusters){
        stop("The number of clusters in 'design' (design$n.clusters) and 'n' (length(n)) do not match.")
      }
      nMat <- matrix(rep(n, each = design$total.time), ncol = design$total.time, byrow = TRUE)#Turns vector into COLUMNS of nMat
    }else if (is.matrix(n)) {
      #n is a matrix with one row per cluster and one column per time point
      if ((nrow(n) != design$n.clusters)|(ncol(n) != design$total.time)){
        stop("The number of clusters and/or time steps in 'design' (design$n.clusters and design$total.time) and 'n' (number of rows and columns, respectively) do not match.")
      }
      nMat <- n#
    }
    Wmat.partial <- matrix(0,nrow=I*J,ncol=I*J)#Wmat.blk doesn't get used anywhere else, so just setting up the Wmat.partial
    for(i in 1:I){#for each of the I total clusters...
      Wmat.partial.i <- tau^2 + diag(sigSq/nMat[i,])+diag(gamma^2,length(nMat[i,]))#Wmat.partial for cluster i
      kronecker.diag.i <- rep(0,I)
      kronecker.diag.i[i] <- 1
      Wmat.partial <- Wmat.partial+kronecker(diag(kronecker.diag.i),Wmat.partial.i)#For each cluster, we fill in another section of the block diagonal matrix
    }
  }
  #Could absolutely do the above without a loop, but efficiency probably not as important here
  Xij.Xil.ARRAY <- array(NA, c(J, J, length(I.rep)))
  for (i.indx in 1:length(I.rep)) {
    Xi.jl <- swDesignUnique[i.indx, ]
    Xij.Xil.ARRAY[, , i.indx] <- (Xi.jl %o% Xi.jl) * eta^2 + outer(Xi.jl, Xi.jl, "+") * rho * eta * tau
  }
  Xij.Xil.LIST_blk <- sapply(1:length(I.rep), function(x) NULL)
  for (i.indx in 1:length(I.rep)) {
    Xij.Xil.LIST_blk[[i.indx]] <- kronecker(diag(1, I.rep[i.indx]), Xij.Xil.ARRAY[, , i.indx])
  }
  W.eta <- blkDiag(Xij.Xil.LIST_blk)
  Wmat <- Wmat.partial + as.matrix(W.eta)
  var.theta.WLS <- solve(t(Xmat) %*% solve(Wmat) %*% Xmat)["X.ij", "X.ij"]
  pwrWLS <- pnorm(abs(theta)/sqrt(var.theta.WLS) - qnorm(1 - alpha/2)) + pnorm(-abs(theta)/sqrt(var.theta.WLS) - qnorm(1 - alpha/2))
  rslt <- pwrWLS
  #Not using the options below for now, but keeping them here in case we want to add them to the shiny app
  if (eta == 0 & gamma == 0 & length(n) == 1) {#we can only calculate closed form when gamma=0 and n is an integer
    X <- swDesign
    U <- sum(X)
    W <- sum(colSums(X)^2)
    V <- sum(rowSums(X)^2)
    sigSq <- sigSq/n
    var.theta.CLOSED <- I * sigSq * (sigSq + J * tau^2)/((I * 
                                                            U - W) * sigSq + (U^2 + I * J * U - J * W - I * V) * 
                                                           tau^2)
    pwrCLOSED <- pnorm(abs(theta)/sqrt(var.theta.CLOSED) - 
                         qnorm(1 - alpha/2)) + pnorm(-abs(theta)/sqrt(var.theta.CLOSED) - 
                                                       qnorm(1 - alpha/2))
  }else {
    pwrCLOSED <- NA
  }
  if (retDATA) 
    rslt <- list(design = design, n = n, mu0 = mu0, mu1 = mu1, 
                 tau = tau, eta = eta, rho = rho, sigma = sigma, gamma=gamma, alpha = alpha, 
                 Xmat = Xmat, Wmat = Wmat, var.theta.WLS = var.theta.WLS, 
                 pwrWLS = pwrWLS, pwrCLOSED = pwrCLOSED)
  rslt
}

####################################################
##User Interface code
####################################################


# Define UI for app
ui <- fluidPage(
  
  # App title ----
  titlePanel("Stepped Wedge Power Calculation"),
  
  #Set up panel format
  tabsetPanel(
  
    ####################################################
    ##Main tab
    ####################################################
    
    #First tab - main arguments and output
    tabPanel("Main results",
             
      # Sidebar for inputs
      sidebarLayout(
    
    
    sidebarPanel(

      #Using wellPanel to get scrolling so that power plot can be seen while each of the sliders is moved
      wellPanel(id = "tPanel",style = "overflow-y:scroll; max-height: 80vh",
                
                
      h5(tags$b("Parameterization"),actionButton("info.variance_param", label=NULL,icon=icon("question-circle"),style='border-radius: 50%; padding: 1px; font-size:80%')),#,width='10%' padding:3px;
      radioButtons(inputId = "variance_param", label=NULL,
                             c("Random effects (more flexible)" = 1, "ICC/CAC (less flexible)" = 2)),      
      

      h5(tags$b("Distribution"),actionButton("info.distn", label=NULL,icon=icon("question-circle"),style='border-radius: 50%; padding: 1px; font-size:80%')),#,width='10%' padding:3px;
      selectInput("distn", label=NULL,#"Distribution",
                  c("Gaussian" = "gaussian",
                    "Binomial" = "binomial")),
      
      h5(tags$b("Number of sequences"),actionButton("info.waves", label=NULL,icon=icon("question-circle"),style='border-radius: 50%; padding: 1px; font-size:80%')),#,width='10%' padding:3px;
      sliderInput(inputId = "waves",
                  label = NULL,#"Number of waves",
                  min = 1,
                  max = 20,
                  value = 5),#Used to be 10
      
      conditionalPanel(#This gets 'outsourced' to another tab if the number of clusters per sequence is not constant
        condition="input.clusters == 1",
      h5(tags$b("Number of clusters per sequence"),actionButton("info.nclust", label=NULL,icon=icon("question-circle"),style='border-radius: 50%; padding: 1px; font-size:80%')),#,width='10%' padding:3px;
      sliderInput(inputId = "nclust",
                  label = NULL,#"Number of clusters per wave",
                  min = 1,
                  max = 20,
                  value = 5)),
      
      
      conditionalPanel(
        condition="input.whichN == 1",
      h5(tags$b("Number of observations per cluster at each time point"),actionButton("info.n", label=NULL,icon=icon("question-circle"),style='border-radius: 50%; padding: 1px; font-size:80%')),#,width='10%' padding:3px;
      sliderInput(inputId = "n",
                  label = NULL,#"Number of observations for each cluster at each time point",
                  min = 1,
                  max = 200,
                  value = 30)),
      

      h5(tags$b("Mean outcome (trt)"),actionButton("info.mu1", label=NULL,icon=icon("question-circle"),style='border-radius: 50%; padding: 1px; font-size:80%')),#,width='10%' padding:3px;
      sliderInput(inputId = "mu1",
                  label = NULL,#"Effect size (trt)",
                  min = -5,#These used to be +/- 10
                  max = 5,
                  value = 0),
      
      h5(tags$b("Mean  outcome (control)"),actionButton("info.mu0", label=NULL,icon=icon("question-circle"),style='border-radius: 50%; padding: 1px; font-size:80%')),#,width='10%' padding:3px;
      sliderInput(inputId = "mu0",
                  label = NULL,#"Effect size (control)",
                  min = -5,
                  max = 5,
                  value = 0),
      
      #Tau, eta, gamma
      conditionalPanel(
        condition="input.variance_param == 1",
        
      h5(tags$b("Tau"),actionButton("info.tau", label=NULL,icon=icon("question-circle"),style='border-radius: 50%; padding: 1px; font-size:80%')),#,width='10%' padding:3px;
      sliderInput(inputId = "tau",
                  label = NULL,#"Tau",
                  min = 0,
                  max = 5,
                  value = 0,step=0.2),
      
      h5(tags$b("Gamma"),actionButton("info.gamma", label=NULL,icon=icon("question-circle"),style='border-radius: 50%; padding: 1px; font-size:80%')),#,width='10%' padding:3px;
      sliderInput(inputId = "gamma",
                  label = NULL,#"Gamma",
                  min = 0,
                  max = 5,
                  value = 0,step=0.2),
      
      h5(tags$b("Eta"),actionButton("info.eta", label=NULL,icon=icon("question-circle"),style='border-radius: 50%; padding: 1px; font-size:80%')),#,width='10%' padding:3px;
      sliderInput(inputId = "eta",
                  label = NULL,#"Eta",
                  min = 0,
                  max = 5,
                  value = 0,step=0.2)),
      
      #ICC and CAC
      
      conditionalPanel(
        condition="input.variance_param == 2",
        
        h5(tags$b("within-period ICC"),actionButton("info.icc", label=NULL,icon=icon("question-circle"),style='border-radius: 50%; padding: 1px; font-size:80%')),#,width='10%' padding:3px;
        sliderInput(inputId = "icc",
                    label = NULL,
                    min = 0,
                    max = 1,
                    value = 0,step=0.1),
        
        h5(tags$b("CAC"),actionButton("info.cac", label=NULL,icon=icon("question-circle"),style='border-radius: 50%; padding: 1px; font-size:80%')),#,width='10%' padding:3px;
        sliderInput(inputId = "cac",
                    label = NULL,
                    min = 0,
                    max = 1,
                    value = 0,step=0.1)
        
        ),
      
      

      
      #Only display sigma option if distn == "gaussian"
      conditionalPanel(
        condition="input.distn == 'gaussian'",
      h5(tags$b("Sigma"),actionButton("info.sigma", label=NULL,icon=icon("question-circle"),style='border-radius: 50%; padding: 1px; font-size:80%')),#,width='10%' padding:3px;
      sliderInput(inputId = "sigma",
                  label = NULL,#"Sigma",
                  min = 0,
                  max = 5,
                  value = 2,step=0.2)),
      
      conditionalPanel(
        condition="input.variance_param == 1",
      h5(tags$b("Rho"),actionButton("info.rho", label=NULL,icon=icon("question-circle"),style='border-radius: 50%; padding: 1px; font-size:80%')),#,width='10%' padding:3px;
      sliderInput(inputId = "rho",
                   label = NULL,#"Rho",
                   min = -1,
                   max = 1,
                   value = 0,step=0.2)),
      
      
      
      h5(tags$b("Alpha"),actionButton("info.alpha", label=NULL,icon=icon("question-circle"),style='border-radius: 50%; padding: 1px; font-size:80%')),#,width='10%' padding:3px;
      numericInput(inputId = "alpha",
                   label = NULL,#"Alpha",
                   min = 0,
                   max = 1,
                   value = 0.05),
      
      
      
      h5(tags$b("Fractional treatment effect"),actionButton("info.tx.effect", label=NULL,icon=icon("question-circle"),style='border-radius: 50%; padding: 1px; font-size:80%')),#,width='10%' padding:3px;
      textInput(inputId = "tx.effect",#There may be a better way to set up this input, but I could not think of a simple way
                label = NULL,"1")
      ),#end Wellpanel
      
      
      
      #turn power curve on and off
      checkboxInput("plotON", 
                    "Plot a power curve",
                    value=FALSE
      ),
      
      #turn study design plot on and off
      checkboxInput("designplotON", 
                    "Visualize study design",
                    value=FALSE
      )
      

      
      
    ),#end of sidebarPanel
    
    # Main panel for displaying outputs
    mainPanel(
      
      verbatimTextOutput("textoutput"),#Power (text output for single value)
      
      busyIndicator(wait = 1500),#Busy indicator for power curve
      plotOutput(outputId = "distPlot"),#power curve
      verbatimTextOutput("textoutput2"),#other outputs describing the design
      plotOutput(outputId = "designPlot")#visualization of design
      

    )  
    )
  ),
  
  ####################################################
  ##Second tab
  ####################################################
  
  #Tab for more complicated, less common  design options
  tabPanel("Additional study design options",h1("Change other elements of the study design")
           
           ,
           
           fluidRow(
             
             
             mainPanel(
               
               h5(tags$b("Special note:")),
               h5("Some of the options on this tab require some additional reactivity with the first tab.  For example, if the number of clusters per sequence varies over time and the number of sequences is changed, the 'error' can be fixed by going back to this tab momentarily."),
               
               
               
               h5(tags$b("Extra time"),actionButton("info.extra.time", label=NULL,icon=icon("question-circle"),style='border-radius: 50%; padding: 1px; font-size:80%')),#,width='10%' padding:3px;
               numericInput(inputId = "extra.time",
                            label = NULL,
                            min = 0,
                            value = 0),
               
               h5(tags$b("Start with control"),actionButton("info.all.ctl.time0", label=NULL,icon=icon("question-circle"),style='border-radius: 50%; padding: 1px; font-size:80%')),#,width='10%' padding:3px;
               selectInput("all.ctl.time0", label=NULL,
                           c("Yes" = "TRUE",
                             "No" = "FALSE")),
               
               
               h5(tags$b("Customize the number of clusters per sequence"),actionButton("info.clusters", label=NULL,icon=icon("question-circle"),style='border-radius: 50%; padding: 1px; font-size:80%')),#,width='10%' padding:3px;
               radioButtons("clusters",label=NULL,c("Constant"=1,"Varying"=2)
               ),
               
               conditionalPanel(
                 condition="input.clusters == 2",
                 h5(tags$b("Number of clusters per sequence"),actionButton("info.clusters2", label=NULL,icon=icon("question-circle"),style='border-radius: 50%; padding: 1px; font-size:80%')),
                 tableOutput("clustersperwave")
               ),
               
               
               h5(tags$b("Customize the number of observations"),actionButton("info.whichN", label=NULL,icon=icon("question-circle"),style='border-radius: 50%; padding: 1px; font-size:80%')),#,width='10%' padding:3px;
               radioButtons("whichN",label=NULL,c("Constant"=1,"Varying by cluster"=2,"Varying by cluster and time"=3)
               ),
               #whichN when number of observations vary over CLUSTER (but not time)
               conditionalPanel(
                 condition="input.whichN == 2",
                 h5(tags$b("Number of observations"),actionButton("info.whichN2", label=NULL,icon=icon("question-circle"),style='border-radius: 50%; padding: 1px; font-size:80%')),
                 tableOutput("Ncluster")
               ),
               #whichN when number of observations is a matrix
               conditionalPanel(
                 condition="input.whichN == 3",
                 h5(tags$b("Number of observations"),actionButton("info.whichN3", label=NULL,icon=icon("question-circle"),style='border-radius: 50%; padding: 1px; font-size:80%')),
                 tableOutput("Ntimecluster")
               )
             )
             
           )
  ),
  
  ####################################################
  ##Third tab
  ####################################################
  
  #Control values in sliders
  tabPanel("Customize sliders",h1("Control range of sliders in first tab"),
           
           fluidRow(
             
             
             sidebarPanel(#Note: it would be possible to hide some of these based on choices in the first tab.
               #If it seemed useful, could expand this to do max, min, and step size for each slider; for now, just including the most essential ones
               

               numericInput(inputId = "waves_max",
                            label = "Maximum number of sequences",
                            min = 1,
                            #max = 20,
                            value = 20),
               numericInput(inputId = "nclust_max",
                           label = "Maximum number of clusters per sequence",
                           min = 1,
                           #max = 20,
                           value = 20),
               numericInput(inputId = "n_max",
                           label = "Maximum number of observations for each cluster at each time point",
                           min = 1,
                           #max = 200,
                           value = 200),
               numericInput(inputId = "mu1_max",
                           label = "Maximum mean outcome (trt)",
                           #min = -10,
                           #max = 10,
                           value = 5),#These used to be +/- 10
               numericInput(inputId = "mu1_min",
                            label = "Minimum mean outcome (trt)",
                            #min = -10,
                            #max = 10,
                            value = -5),
               numericInput(inputId = "mu0_max",
                           label = "Maximum mean outcome (control)",
                           #min = -10,
                           #max = 10,
                           value = 5),
               numericInput(inputId = "mu0_min",
                            label = "Minimum mean outcome (control)",
                            #min = -10,
                            #max = 10,
                            value = -5),
               numericInput(inputId = "tau_max",
                           label = "Maximum tau",
                           min = 0,
                           #max = 5,
                           value = 5),
               numericInput(inputId = "tau_step",
                            label = "Size of tau intervals",
                            min = 0,
                            #max = 5,
                            value = 0.2),
               numericInput(inputId = "gamma_max",
                            label = "Maximum gamma",
                            min = 0,
                            #max = 5,
                            value = 5),
               numericInput(inputId = "gamma_step",
                            label = "Size of gamma intervals",
                            min = 0,
                            #max = 5,
                            value = 0.2),
               numericInput(inputId = "eta_max",
                           label = "Maximum eta",
                           min = 0,
                           #max = 5,
                           value = 5),
               numericInput(inputId = "eta_step",
                            label = "Size of eta intervals",
                            min = 0,
                            #max = 5,
                            value = 0.2),
               numericInput(inputId = "sigma_max",
                           label = "Maximum sigma",
                           min = 0,
                           #max = 5,
                           value = 5),
               numericInput(inputId = "sigma_step",
                            label = "Size of sigma intervals",
                            min = 0,
                            #max = 5,
                            value = 0.2),
               numericInput(inputId = "rho_max",
                            label = "Maximum rho",
                            min = -1,
                            max = 1,
                            value = 1),
               numericInput(inputId = "rho_min",
                            label = "Minimum rho",
                            min = -1,
                            max = 1,
                            value = -1),
               numericInput(inputId = "rho_step",
                            label = "Size of rho intervals",
                            min = 0,
                            max = 2,
                            value = 0.2),
               numericInput(inputId = "icc_max",
                            label = "Maximum ICC",
                            min = 0,
                            max = 1,
                            value = 1),
               numericInput(inputId = "icc_min",
                            label = "Minimum ICC",
                            min = 0,
                            max = 1,
                            value = 0),
               numericInput(inputId = "cac_max",
                            label = "Maximum CAC",
                            min = 0,
                            max = 1,
                            value = 1),
               numericInput(inputId = "cac_min",
                            label = "Minimum CAC",
                            min = 0,
                            max = 1,
                            value = 0)
             )
             
           )
  ),
  
  ####################################################
  ##Fourth tab
  ####################################################
  
  #Tab to customize power curve
  tabPanel("Customize plots",h1("Control characteristics of the plots in first tab"),
           
           fluidRow(
             
             
             mainPanel(
               
               h3("Power curve"),
               
               
               h5(tags$b("Strategy for defining width of possible effect sizes"),actionButton("info.effectsize_mode", label=NULL,icon=icon("question-circle"),style='border-radius: 50%; padding: 1px; font-size:80%')),#,width='10%' padding:3px;
               selectInput("effectsize_mode", label=NULL,
                           c("Fixed maximum and minimum" = 1,
                             "Range around effect size (trt-control)" = 2)),
               #If using a fixed max and min...
               conditionalPanel(
                 condition="input.effectsize_mode == 1",
               numericInput(inputId = "effectsize_min",
                            label = "Minimum effect size",
                            value = -2),
               numericInput(inputId = "effectsize_max",
                            label = "Maximum effect size",
                            value = 2)),
               #If using a range...
               conditionalPanel(
                 condition="input.effectsize_mode == 2",
               numericInput(inputId = "effectsize_range",
                            label = "Effect size range",
                            min = 0,
                            value = 2)),#Used to be 3
               
               numericInput(inputId = "effectsize_distance",
                            label = "Distance in effect size between plotted points",
                            min = 0,
                            value = 0.1),
               
               selectInput("plot_style", label="Style of plot",
                           c("Points only" = "points",
                             "Lines only" = "lines",
                             "Both points and lines" = "both")),
               
               #Download points on power curve (May be a problem with reactivity here if user does not return to the first tab after making changes)
               downloadButton("downloadData", "Download"),
               
               h3("Study design plot"),
               selectInput("plot_extra_lines", label="Scale of design",
                           c("Subdivide sequences by cluster" = "cluster",
                             "Plot only sequences" = "sequence"))
             )
             
           )
  ),

  ####################################################
  ##Fifth tab
  ####################################################
  
  tabPanel("Help and References",
           
           fluidRow(
             
             
             mainPanel(
               
               h4("Report bugs not addressed below to: voldal@uw.edu"),
               tags$p("Contributors (alphabetical order): Navneet Hakhu, Patrick Heagerty, Jim Hughes, Emily Voldal"),
               tags$p("Last updated: 8/28/2019"),
               tags$p("The full swCRTdesign R package can be found on CRAN (https://cran.r-project.org/web/packages/swCRTdesign/index.html), and the code for this Shiny app is available on Github (https://github.com/swCRTdesign/Stepped-wedge-power-calculation)."),
               tags$p("Research reported in this work was partially funded through a Patient-Centered Outcomes Research Institute (PCORI) Award (ME-1507-31750).  The statements in this work are solely the responsibility of the authors and do not necessarily represent the views of the Patient-Centered Outcomes Research Institute (PCORI), its Board of Governors, or Methodology Committee."),
               
               
               includeHTML("manual_markdown.html")#Help manual
               #If for some reason you are having issues running the app with the manual (or you don't have this manual file), just comment out the line above.
             )
             
           )
  )

))

####################################################
##Server code
####################################################



server <- function(input, output,session) {
  
  ####################################################
  ##Update sliders in main tab based on user inputs
  ####################################################
  
  
  observe({
    
    
    updateSliderInput(session, "waves",
                      max = input$waves_max)
    updateSliderInput(session, "nclust",
                      max = input$nclust_max)
    updateSliderInput(session, "n",
                      max = input$n_max)
    updateSliderInput(session, "mu1",
                      min = input$mu1_min,
                      max = input$mu1_max,
                      step=(input$mu1_max-input$mu1_min)/20)
    updateSliderInput(session, "mu0",
                      min = input$mu0_min,
                      max = input$mu0_max,
                      step=(input$mu0_max-input$mu0_min)/20)
    updateSliderInput(session, "tau",
                      max = input$tau_max,step = input$tau_step)
    updateSliderInput(session, "eta",
                      max = input$eta_max,step = input$eta_step)
    updateSliderInput(session, "gamma",
                      max = input$gamma_max,step = input$gamma_step)
    updateSliderInput(session, "sigma",
                      max = input$sigma_max,step = input$sigma_step)
    updateSliderInput(session, "rho",
                      min = input$rho_min,
                      max = input$rho_max,
                      step= input$rho_step)
    updateSliderInput(session, "icc",
                      min = input$icc_min,
                      max = input$icc_max,
                      step=(input$icc_max-input$icc_min)/20)
    updateSliderInput(session, "cac",
                      min = input$cac_min,
                      max = input$cac_max,
                      step=(input$cac_max-input$cac_min)/20)
    
    

  })
  
  observeEvent({input$distn == "binomial"},{
    #if distribution is binomial, average effects can only be in (0,1)
    #But still give user ability to adjust by only changing this once when dist'n is changed
    if (input$distn == "binomial" ){
      updateNumericInput(session, "mu1_min",
                         value=0)
      updateNumericInput(session, "mu1_max",
                         value=1)
      updateNumericInput(session, "mu0_min",
                         value=0)
      updateNumericInput(session, "mu0_max",
                         value=1)
    }
  })
  
  ####################################################
  ##Help button text
  ####################################################
  
  
  #Help panels
  observeEvent(input$info.waves, {
    showModal(modalDialog(
      title = "Number of sequences",
      "All clusters that start the intervention at a given time point are collectively referred to as a sequence.  For simple designs, the number of sequences will be the number of time points where data is collected minus one (since all clusters start in the control group).",
      easyClose = TRUE
    ))
  })
  observeEvent(input$info.nclust, {
    showModal(modalDialog(
      title = "Number of clusters per sequence",
      "The number of clusters that switch to treatment at each sequence.  All clusters that start the intervention at a given time point are collectively referred to as a sequence.  For designs where the number of clusters per sequence varies, see other tabs.",
      easyClose = TRUE
    ))
  })
  observeEvent(input$info.n, {
    showModal(modalDialog(
      title = "Number of observations per cluster at each time point",
      "The number of observations for each cluster at each time point (the same for all clusters at all time points).  For designs where the number of observations varies, see other tabs.",
      easyClose = TRUE
    ))
  })
  observeEvent(input$info.mu1, {
    showModal(modalDialog(
      title = "Mean outcome (trt)",
      "Mean of the outcome in the treatment group.",
      easyClose = TRUE
    ))
  })
  observeEvent(input$info.mu0, {
    showModal(modalDialog(
      title = "Mean outcome (control)",
      "Mean of the outcome in the control group.  Note that the treatment effect is the difference between the treatment average and the control average.",
      easyClose = TRUE
    ))
  })
  observeEvent(input$info.tau, {
    showModal(modalDialog(
      title = "Tau",
      "The standard deviation of the random intercepts.",
      easyClose = TRUE
    ))
  })
  observeEvent(input$info.eta, {
    showModal(modalDialog(
      title = "Eta",
      "The standard deviation of the random treatment effects.  Setting eta to 0 corresponds to a model with a fixed treatment effect.",
      easyClose = TRUE
    ))
  })
  observeEvent(input$info.gamma, {
    showModal(modalDialog(
      title = "Gamma",
      "The standard deviation of the random time effects.  Setting gamma to 0 corresponds to a model with a fixed time effect.",
      easyClose = TRUE
    ))
  })
  observeEvent(input$info.sigma, {
    showModal(modalDialog(
      title = "Sigma",
      "The standard deviation of the random error (for a Gaussian distribution only).",
      easyClose = TRUE
    ))
  })
  observeEvent(input$info.icc, {
    showModal(modalDialog(
      title = "within-period ICC",
      "The within-period intra-cluster correlation.  This is the correlation between individuals sampled from the same cluster at the same time.  Equivalently, this is the ratio of between-cluster variance to total variance.",
      easyClose = TRUE
    ))
  })
  observeEvent(input$info.cac, {
    showModal(modalDialog(
      title = "CAC",
      "The cluster auto-correlation.  This is the correlation between two population means from the same cluster at different times.  This is calculated as the ratio of between-period ICC to within-period ICC.",
      easyClose = TRUE
    ))
  })
  observeEvent(input$info.distn, {
    showModal(modalDialog(
      title = "Distribution",
      "Assumed distribution of the outcome.  Choose 'Gaussian' for continuous data, and 'binomial' for binary data.  Note that for the binomial distribution, the parameter 'sigma' is determined by the treatment effect.",
      easyClose = TRUE
    ))
  })
  observeEvent(input$info.alpha, {
    showModal(modalDialog(
      title = "Alpha",
      "Two-sided statistical significance level.  Default is 0.05.",
      easyClose = TRUE
    ))
  })
  observeEvent(input$info.rho, {
    showModal(modalDialog(
      title = "Rho",
      "Correlation between random intercepts and random treatment effects.  Default is 0 (no correlation).",
      easyClose = TRUE
    ))
  })
  observeEvent(input$info.tx.effect, {
    showModal(modalDialog(
      title = "Fractional treatment effect",
      "Fractional treatment effect upon crossing over from control.  A value of 1 corresponds to the standard treatment effect.  Enter a list of numbers (separated by commas) between 0 and 1.  Each number in the list corresponds to the fractional treatment effect upon crossing over from control for a different time point.
      If the list is shorter than the total number of time points after crossing over, the remaining time points will have a fractional treatment effect of 1.  If the vector is longer than the total number of time points after crossing over, not all elements in the vector will be used.  The default value is 1.  For example, if the fractional treatment effect is '0.5,0.75', then the first time a cluster receives treatment the effect will be half the size of the specified treatment effect, the next time point it will have 75% of the specified treatment effect, and after that the cluster will experience the full treatment effect. ",
      easyClose = TRUE
    ))
  })
  observeEvent(input$info.variance_param, {
    showModal(modalDialog(
      title = "Parameterization",
      "How to parameterize the variance in the data, which may be correlated in several different ways; see help tab.  Note that using the ICC/CAC parameterization requires assuming that there are no random treatment effects.",
      easyClose = TRUE
    ))
  })
  observeEvent(input$info.extra.time, {
    showModal(modalDialog(
      title = "Extra time",
      "Number of additional time steps beyond the standard design - that is, additional observations after all clusters have been switched to treatment.",
      easyClose = TRUE
    ))
  })
  observeEvent(input$info.all.ctl.time0, {
    showModal(modalDialog(
      title = "Start with control",
      "All clusters receive control at the first time point.  If 'no', clusters in the first sequence receive the treatment at the first time point, and never the control.",
      easyClose = TRUE
    ))
  })
  observeEvent(input$info.clusters, {
    showModal(modalDialog(
      title = "Customize the number of clusters per sequence",
      "The default ('constant') is that the same number of clusters switch to treatment in each sequence.  Alternatively ('varying'), the number of clusters that switch to treatment in each sequence can vary over time.",
      easyClose = TRUE
    ))
  })
  observeEvent(input$info.whichN, {
    showModal(modalDialog(
      title = "Customize the number of observations",
      "The default ('constant') is that there is the same number of observations in each cluster at each time point.  Alternatively ('varying by cluster'), the number of observations can be different for each cluster (but be the same at every time point within each cluster).  Or ('varying by cluster and time') the number of observations can be different for each cluster at each time point.",
      easyClose = TRUE
    ))
  })
  observeEvent(input$info.clusters2, {
    showModal(modalDialog(
      title = "Number of clusters per sequence",
      "Identify how many clusters change from control to intervention at each sequence.  A value of 0 means that no clusters introduce the intervention at a given time.",
      easyClose = TRUE
    ))
  })
  observeEvent(input$info.whichN2, {
    showModal(modalDialog(
      title = "Number of observations",
      "Number of observations for each cluster at every time point.  That is, the number of observations for a specific cluster is constant over time.  Each row represents one cluster, and rows are ordered by the order in which the clusters receive treatment.",
      easyClose = TRUE
    ))
  })
  observeEvent(input$info.whichN3, {
    showModal(modalDialog(
      title = "Number of observations",
      "Number of observations for each cluster at each time point.  The columns correspond to the times, and the rows correspond to the clusters.  For example, the number in the first column and third row would be the number of observations in the third cluster at the first time point.  The order of the rows should correspond to the order in which the clusters switch to treatment, with the earlier sequences in the upper rows.",
      easyClose = TRUE
    ))
  })
  #info.effectsize_mode
  observeEvent(input$info.effectsize_mode, {
    showModal(modalDialog(
      title = "Strategy for defining width of possible effect sizes",
      "This is the method which will determine the range of the x-axis in the power plot.  If 'Range around effect size' is selected, the effect size defined in the main panel will determine the center of the plot (indicated by a vertical line).  For a binomial distribution, the mean outcome in the control group on the 'Main results' tab is used to generate the plot, and changes in effect size correspond to changes in the treatment group mean outcome.",
      easyClose = TRUE
    ))
  })
  
  ####################################################
  ##Reactive tables
  ####################################################
  
  
  #Reactive table for the number of clusters per wave
  output$clustersperwave <-renderTable({
    num.inputs.col1 <- paste0("<input id='waves.r", 1, "c", 1, "' class='shiny-bound-input' type='number' value='1'>")
    df <- data.frame(num.inputs.col1)
    if (input$waves >= 2){#want input$waves
      for (i in 2:input$waves){
        num.inputs.coli <- paste0("<input id='waves.r", 1, "c", i, "' class='shiny-bound-input' type='number' value='1'>")
        df <- cbind(df,num.inputs.coli)
      }
    }
    
    colnames(df) <- paste0("Sequence ",as.numeric(1:input$waves))
    df
  }, sanitize.text.function = function(x) x)
  
  #Reactive table for the number of observations (time- and cluster-varying)
  output$Ntimecluster <-renderTable({
    #calculate total time: number of waves + extra.time + all.ctl.time0
    ncol2 <- input$waves+input$extra.time+as.numeric(as.logical(input$all.ctl.time0))
    #calculate number of clusters: depends on 'clusters'
    if (input$clusters == 1){
      nrow2 <- input$nclust*input$waves#number of clusters = number of clusters per wave * number of waves
    }
    if (input$clusters == 2){
      #Need to sum up entries in the vector for number of clusters
      nrow2 <- 0
      for (i in 1:input$waves){
        eval(parse(text=paste0("nrow2 <- nrow2 + input$waves.r1c",i)))
      }
    }
    
    num.inputs.col1 <- paste0("<input id='Ntimecluster.r", 1:nrow2, "c", 1, "' class='shiny-bound-input' type='number' value='1'>")
    df <- data.frame(num.inputs.col1)
    if (ncol2 >= 2){
      for (i in 2:ncol2){
        num.inputs.coli <- paste0("<input id='Ntimecluster.r", 1:nrow2, "c", i, "' class='shiny-bound-input' type='number' value='1'>")
        df <- cbind(df,num.inputs.coli)
      }
    }
    
    colnames(df) <- paste0("Time ",as.numeric(1:ncol2))
    rownames(df) <- paste0( 1:nrow2)
    df
  }, sanitize.text.function = function(x) x,rownames=TRUE)
  
  #Reactive table for the number of observations (only time-varying) - NOT USING THIS, but keeping it here because it would be a valid option
  # output$Ntime <-renderTable({
  #   #calculate total time: number of waves + extra.time + all.ctl.time0
  #   ncol3 <- input$waves+input$extra.time+as.numeric(as.logical(input$all.ctl.time0))
  #   
  #   num.inputs.col1 <- paste0("<input id='Ntime.r", 1, "c", 1, "' class='shiny-bound-input' type='number' value='1'>")
  #   df <- data.frame(num.inputs.col1)
  #   if (ncol3 >= 2){#want input$waves
  #     for (i in 2:ncol3){
  #       num.inputs.coli <- paste0("<input id='Ntime.r", 1, "c", i, "' class='shiny-bound-input' type='number' value='1'>")
  #       df <- cbind(df,num.inputs.coli)
  #     }
  #   }
  #   
  #   colnames(df) <- paste0("Time ",as.numeric(1:ncol3))
  #   df
  # }, sanitize.text.function = function(x) x)
  
  #Reactive table for the number of observations per cluster (not time varying)
  output$Ncluster <-renderTable({
    #calculate number of clusters: depends on 'clusters'
    if (input$clusters == 1){
      nrow4 <- input$nclust*input$waves#number of clusters = number of clusters per wave * number of waves
    }
    if (input$clusters == 2){
      #Need to sum up entries in the vector for number of clusters
      nrow4 <- 0
      for (i in 1:input$waves){
        eval(parse(text=paste0("nrow4 <- nrow4 + input$waves.r1c",i)))
      }
    }
    
    num.inputs.col1 <- paste0("<input id='Ncluster.r", 1:nrow4, "c", 1, "' class='shiny-bound-input' type='number' value='1'>")
    df <- data.frame(num.inputs.col1)
    #if (ncol2 >= 2){#want input$waves
    #  for (i in 2:ncol2){
    #    num.inputs.coli <- paste0("<input id='Ncluster.r", 1:nrow2, "c", i, "' class='shiny-bound-input' type='number' value='1'>")
    #    df <- cbind(df,num.inputs.coli)
    #  }
    #}
    
    #colnames(df) <- paste0("Time ",as.numeric(1:ncol2))
    colnames(df) <- paste0("Number of observations per time point")
    rownames(df) <- paste0( 1:nrow4)
    df
  }, sanitize.text.function = function(x) x,rownames=TRUE)
  
  ####################################################
  ##Plot study design
  ####################################################
  
  
  #Plot of study design
  output$designPlot <- renderPlot({
    if(input$designplotON == TRUE){
      try({
      #Need swdsn object
      tx.effect.num <- as.numeric(unlist(strsplit(input$tx.effect,",")))
      design.vector <- rep(1,length=input$waves)*input$nclust
      #Design vector for the case where the number of clusters varies by wave:
      if (input$clusters == 2){
        for (i in 1:input$waves){
          eval(parse(text=paste0("design.vector[i] <- input$waves.r1c",i)))
        }
      }
      
      design <- swDsn(clusters=design.vector,  extra.time=input$extra.time, all.ctl.time0=as.logical(input$all.ctl.time0), tx.effect=tx.effect.num)
      
      #Plot with ggplot2
      matrix <- design$swDsn.unique.clusters#Design matrix straight from swCRTdesign
      ntime <- ncol(matrix)
      nwave <- nrow(matrix)
      
      matrixlong <- matrix(NA,nrow=nwave*ntime,ncol=3)#Var1,Var2,value
      matrixlong[,1] <- rep(1:nwave,ntime)#waves
      matrixlong[,2] <- rep(1:ntime,each=nwave)#times
      
      #Loop over rows of matrixlong
      for (i in 1:nrow(matrixlong)){
        matrixlong[i,3] <- matrix[matrixlong[i,1],matrixlong[i,2]]
      }
      
      matrixdf <- data.frame(matrixlong)
      names(matrixdf) <- c("Var1","Var2","value")
      matrixdf$value <- as.factor(matrixdf$value)
      
      #make vectors for grid
      timevector <- c(1:ntime,(ntime+1))-0.5
      wavevector <- c(1:nwave,(nwave+1))-0.5
      
      clustervector <- c()
      design.vector.no0 <- design.vector[design.vector != 0]
      for (i in 1:length(design.vector.no0)){
        clustervector <- c(clustervector,seq(from=wavevector[i],to=wavevector[i+1],by=1/(design.vector.no0[i])))
      }
      
      design.plot.base <- ggplot(matrixdf, aes(x = Var2, y = Var1)) + 
        geom_raster(aes(fill=value)) + 
        scale_fill_manual(values=c("white", "grey"),name="Status",labels=c("Control","Treatment")) +labs(x="Time",y="Sequence")+scale_y_reverse(expand=c(0,0),breaks=(1:nwave))+theme_bw()+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.box.background = element_rect(fill = NA,color=NA),legend.key = element_rect(fill = NA, colour = "black",linetype="solid",size=1),plot.title=element_text(hjust=0.5))+
        ggtitle("Study Design")+geom_vline(xintercept=timevector)+geom_hline(yintercept=wavevector)+
        scale_x_continuous(expand=c(0,0),breaks=(1:ntime))
      

      if (input$plot_extra_lines == "cluster"){
        design.plot.base+geom_hline(yintercept=clustervector,linetype="dotted")
      }else{design.plot.base}
      
      
      }, silent=TRUE)
      
      #tryCatch(, error=function(e) print("The study design plot could not be generated.  Double-check inputs, and consult the manual for extra tips."))
    }
  })
  
  ####################################################
  ##Translate from ICC/CAC to random effects
  ####################################################
  
  
  #Take input in either linear model or ICC/CAC form, and return tau, eta, gamma, rho (sigma stays the same) (also depends on binomial/gaussian)
  variance_components_vector <- reactive({
    if (input$variance_param == 1){
      vector <- c(input$tau,input$gamma,input$eta,input$rho)
    }
    if (input$variance_param == 2){
      #Assume eta=0 and rho=0
      if (input$distn == 'gaussian'){
        sigmasq.temp <- input$sigma^2
      }
      if (input$distn != 'gaussian'){
        mubar.temp <- (input$mu1+input$mu0)/2
        sigmasq.temp <- mubar.temp*(1-mubar.temp)
      }
      
      if(input$cac == 1){
          gamma.temp <- 0
          tau.temp <- sqrt(sigmasq.temp*input$icc/(1-input$icc))
      }else{
          gamma.temp <- sqrt(input$icc*sigmasq.temp*(1-input$cac)/(1-input$icc))
          tau.temp <- sqrt(gamma.temp^2*input$cac/(1-input$cac))
        }
      
      vector <- c(tau.temp,gamma.temp,0,0)
    }
    return(vector)#tau,gamma,eta,rho
  })
  
  ####################################################
  ##Power plot data
  ####################################################
  
  
  
  #Set up points from plot to download, and then use them in the plot itself.
  data_for_download <- reactive({ 
    if(input$plotON == TRUE){
    try({
      
      vector.temp <- variance_components_vector()
      input.tau <- vector.temp[1]
      input.gamma <- vector.temp[2]
      input.eta <- vector.temp[3]
      input.rho <- vector.temp[4]
      
      #Process input from the tx.effect list
      tx.effect.num <- as.numeric(unlist(strsplit(input$tx.effect,",")))
      
      design.vector <- rep(1,length=input$waves)*input$nclust
      #Design vector for the case where the number of clusters varies by wave:
      if (input$clusters == 2){
        for (i in 1:input$waves){
          eval(parse(text=paste0("design.vector[i] <- input$waves.r1c",i)))
        }
      }
      
      design <- swDsn(clusters=design.vector,  extra.time=input$extra.time, all.ctl.time0=as.logical(input$all.ctl.time0), tx.effect=tx.effect.num)
      
      #Take input about the number of observations, and return a scalar, vector, or matrix.
      if (input$whichN == 1){#scalar
        flexn <- input$n
      } else if(input$whichN == 2){#vector (varies by cluster)
        #Identify dimension of vector
        #calculate number of clusters: depends on 'clusters'
        if (input$clusters == 1){
          vectorlength <- input$nclust*input$waves#number of clusters = number of clusters per wave * number of waves
        }
        if (input$clusters == 2){
          #Need to sum up entries in the vector for number of clusters
          vectorlength <- 0
          for (i in 1:input$waves){
            eval(parse(text=paste0("vectorlength <- vectorlength + input$waves.r1c",i)))
          }
        }
        #build empty vector
        flexn <- rep(NA,vectorlength)
        #fill it with HTML cells
        for (i in 1:vectorlength){
          eval(parse(text=paste0("flexn[", i, "] <- input$Ncluster.r",i,"c1")))
        }
      } else if(input$whichN == 3){#matrix
        #Identify dimension of matrix
        #calculate total time: number of waves + extra.time + all.ctl.time0
        matrixncol <- input$waves+input$extra.time+as.numeric(as.logical(input$all.ctl.time0))
        #calculate number of clusters: depends on 'clusters'
        if (input$clusters == 1){
          matrixnrow <- input$nclust*input$waves#number of clusters = number of clusters per wave * number of waves
        }
        if (input$clusters == 2){
          #Need to sum up entries in the vector for number of clusters
          matrixnrow <- 0
          for (i in 1:input$waves){
            eval(parse(text=paste0("matrixnrow <- matrixnrow + input$waves.r1c",i)))
          }
        }
        #build empty matrix
        flexn <- matrix(NA,ncol=matrixncol,nrow=matrixnrow)
        #fill it in with HTML cells
        for (i in 1:matrixncol){#For each column i...
          for (j in 1:matrixnrow){#For each row j...
            eval(parse(text=paste0("flexn[", j, ",",i,"] <- input$Ntimecluster.r",j,"c",i)))
          }
        }
      }
      
      GetPower <- function(effectsize){#,input.tau=input.tau,input.eta=input.eta,input.gamma=input.gamma,input.rho=input.rho
        
        if(input$distn == "gaussian"){
          obj <- swPwr.FlexN(design=design,distn=input$distn,n=flexn,mu0=0,mu1=effectsize,tau=input.tau,eta=input.eta,sigma=input$sigma,gamma=input.gamma,rho=input.rho,alpha=input$alpha)
        
          }
        
        if(input$distn == "binomial"){#mu0 is defined by user, and mu1 is mu0+effectsize
          obj <- swPwr.FlexN(design=design,distn=input$distn,n=flexn,mu0=input$mu0,mu1=(input$mu0+effectsize),tau=input.tau,eta=input.eta,gamma=input.gamma,rho=input.rho,alpha=input$alpha,sigma=input$sigma)#removed ,sigma=input$sigma
          
        }
        
        return(obj)
      }
      
      effectsize <- input$mu1-input$mu0
      #Can choose x values with a range or min/max
      if (input$effectsize_mode == 2){
        effectsizelist <- seq(effectsize-input$effectsize_range,effectsize+input$effectsize_range,by=input$effectsize_distance)
      }else{
        effectsizelist <- seq(input$effectsize_min,input$effectsize_max,by=input$effectsize_distance)
      }
      
      if(input$distn == "binomial"){
        #check all the restrictions for mu0 and mu1 (mu0 should be checked elsewhere)
        good_mu1 <- (input$mu0 + effectsizelist) <= 1 & (input$mu0 + effectsizelist) >= 0
        good_var <- input.tau^2+input.eta^2+input.gamma^2 < (effectsizelist+2*input$mu0)/2*(1-(effectsizelist+2*input$mu0)/2)
        effectsizelist <- effectsizelist[good_mu1 & good_var]
      }
      
      

      powerlist <- sapply(effectsizelist,GetPower)
      
      
      data.frame(effectsizelist,powerlist)
      
      
    },silent=TRUE)
    }
    })
  
  ####################################################
  ##Make power plot
  ####################################################
  
  
  output$distPlot <- renderPlot({
    
    
    data_for_plots <- data_for_download()
    #Function to take an effect size and calculate the power given these specific inputs
    #(design is same for all effect sizes)
    
    if(input$plotON == TRUE){
      validate(need(nrow(data_for_plots) > 1, "To plot a power curve with more than one point, increase the effect size range or decrease the distance between points (see 'Customize plots' tab)."))
      
    try({
      
    effectsize <- input$mu1-input$mu0
    
    if (input$plot_style == "points"){
      plot(data_for_plots$effectsizelist,data_for_plots$powerlist,xlab="Effect size",ylab="Power",main="Power Curve")
      abline(v=effectsize)
    }
    if (input$plot_style == "lines"){
      plot(data_for_plots$effectsizelist,data_for_plots$powerlist,xlab="Effect size",ylab="Power",main="Power Curve",type="l")
      abline(v=effectsize)
    }
    if (input$plot_style == "both"){
      plot(data_for_plots$effectsizelist,data_for_plots$powerlist,xlab="Effect size",ylab="Power",main="Power Curve",pch=16)
      abline(v=effectsize)
      lines(data_for_plots$effectsizelist,data_for_plots$powerlist)
    }
    
    
    
    },silent=TRUE)
    }
    if(input$plotON == FALSE){
      plot(0,type='n',main="Select options to the left to plot a power curve.",axes=FALSE,xlab="",ylab="")
    }
    
  })
  

  ####################################################
  ##Make CSV for downloading data
  ####################################################
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("power_plot_points", ".csv", sep = "")
    },#"power_plot_points.csv",
    content = function(file) {
      write.csv(data_for_download(), file, row.names = FALSE)
    }
  )
  
  ####################################################
  ##Calculate power for text output
  ####################################################
  
  output$textoutput <-  renderText({
    
    vector.temp <- variance_components_vector()
    input.tau <- vector.temp[1]
    input.gamma <- vector.temp[2]
    input.eta <- vector.temp[3]
    input.rho <- vector.temp[4]
    
    validate(
      #Main results
      need(input$sigma != 0 | input$variance_param == 1, "When sigma is 0, the random effects parameterization must be used."),
      need((input$sigma != 0) | (input.gamma != 0) | input$distn == "binomial", "The selected variance components correspond to deterministic data; if you are using a Gaussian distribution, try making sigma nonzero."),
      need(input$distn != "binomial" | input.tau^2+input.eta^2+input.gamma^2 < (input$mu1+input$mu0)/2*(1-(input$mu1+input$mu0)/2),"For a binomial outcome, the sum tau^2 + eta^2 + gamma^2 needs to be less than mu*(1-mu), where mu is the average of the mean outcomes from the treatment group and the control group, (mu_0+mu_1)/2."),
      need(input$waves > 1, "The number of sequences must be at least two."),
      need(input$alpha > 0 & input$alpha < 1, "Alpha must be a numeral strictly between 0 and 1."),
      need(all(as.numeric(unlist(strsplit(input$tx.effect,","))) >= 0) & all(as.numeric(unlist(strsplit(input$tx.effect,","))) <= 1), "The fractional treatment effect must be a numeral, or list of numerals separated by commas.  All the values must be between 0 and 1."),
      need(input.rho == 0 |(input.tau != 0 & input.eta != 0),"If either tau or eta is zero, rho must be zero, since fixed effects cannot be correlated.  "),
      need(input$icc != 1 | input$variance_param == 1,"There are multiple combinations of random effects that can make the ICC be 1; if you believe this is a realistic scenario, use the random effect parameterization."),
      
      #Additional study design
      need(input$extra.time >= 0 & input$extra.time%%1 == 0, "Extra time must be a non-negative integer."),
      #I am not checking the two matrix inputs here - somewhat difficult to check, and seems like a pretty low risk of user-input issues.
      
      #Customize sliders
      need(input$waves_max > 1 & input$waves_max%%1 == 0, "The maximum number of sequences must be an integer greater than 1."),
      need(input$nclust_max > 0 & input$nclust_max%%1 == 0, "The maximum number of clusters per sequence must be a positive integer."),
      need(input$n_max > 0 & input$n_max%%1 == 0, "The maximum number of observations must be a positive integer."),
      need(input$mu1_min <= input$mu1_max & input$mu0_min <= input$mu0_max,"The maximum and minimum mean outcomes in both the treatment and control groups need to be numerals.  In addition, each minimum must be at least as small as the corresponding maximum."),
      need(input$distn != "binomial" | (input$mu1_max <= 1 & input$mu0_max <= 1 & input$mu1_min >= 0 & input$mu0_min >= 0),"For a binomial outcome, the mean outcomes must be between 0 and 1."),
      need(input$tau_max >= 0 & input$eta_max >= 0 & input$gamma_max >= 0 & input$sigma_max >= 0, "The maximum values of tau, eta, gamma, and sigma must all be non-negative numerals."),
      need(input$tau_step >= 0 & input$eta_step >= 0 & input$gamma_step >= 0 & input$sigma_step >= 0 & input$rho_step >= 0, "The sizes of the intervals in the sliders for tau, eta, gamma, sigma, and rho must all be non-negative numerals."),
      need(input$rho_min <= input$rho_max & -1 <= input$rho_min & input$rho_max <= 1,"The maximum and minimum values of rho need to be numerals between -1 and 1.  In addition, the minimum must be at least as small as the maximum."),
      need(0 <= input$icc_min & input$icc_min <= input$icc_max & input$icc_max <= 1, "The maximum and minimum values of the ICC need to be numerals between 0 and 1.  In addition, the minimum must be at least as small as the maximum."),
      need(0 <= input$cac_min & input$cac_min <= input$cac_max & input$cac_max <= 1, "The maximum and minimum values of the CAC need to be numerals between 0 and 1.  In addition, the minimum must be at least as small as the maximum."),
      
      #Customize curve
      need((is.na(input$effectsize_min) & is.na(input$effectsize_max)) | input$effectsize_min <= input$effectsize_max ,"When customizing the power curve, the minimum effect size needs to be at least as small as the maximum effect size."),
      need(input$effectsize_range >= 0, "When customizing the power curve, the range of possible effect sizes needs to be a non-negative numeral."),
      need(input$effectsize_distance > 0, "When customizing the power curve, the distance between points needs to be a positive numeral.")
    )
    

    
    tx.effect.num <- as.numeric(unlist(strsplit(input$tx.effect,",")))
    design.vector <- rep(1,length=input$waves)*input$nclust#This is the design vector for the simple case
    #Design vector for the case where the number of clusters varies by wave:
    if (input$clusters == 2){
    for (i in 1:input$waves){
      eval(parse(text=paste0("design.vector[i] <- input$waves.r1c",i)))
    }
    }
    
    design <- swDsn(clusters=design.vector,  extra.time=input$extra.time, all.ctl.time0=as.logical(input$all.ctl.time0), tx.effect=tx.effect.num)
    
    #Assign n as a scalar, vector, or matrix
    #Take input about the number of observations, and return a scalar, vector, or matrix.
    if (input$whichN == 1){#scalar
      flexn <- input$n
    } else if(input$whichN == 2){#vector (varies by cluster)
      #Identify dimension of vector
      #calculate number of clusters: depends on 'clusters'
      if (input$clusters == 1){
        vectorlength <- input$nclust*input$waves#number of clusters = number of clusters per wave * number of waves
      }
      if (input$clusters == 2){
        #Need to sum up entries in the vector for number of clusters
        vectorlength <- 0
        for (i in 1:input$waves){
          eval(parse(text=paste0("vectorlength <- vectorlength + input$waves.r1c",i)))
        }
      }
      #build empty vector
      flexn <- rep(NA,vectorlength)
      #fill it with HTML cells
      for (i in 1:vectorlength){
        eval(parse(text=paste0("flexn[", i, "] <- input$Ncluster.r",i,"c1")))
      }
    } else if(input$whichN == 3){#matrix
      #Identify dimension of matrix
      #calculate total time: number of waves + extra.time + all.ctl.time0
      matrixncol <- input$waves+input$extra.time+as.numeric(as.logical(input$all.ctl.time0))
      #calculate number of clusters: depends on 'clusters'
      if (input$clusters == 1){
        matrixnrow <- input$nclust*input$waves#number of clusters = number of clusters per wave * number of waves
      }
      if (input$clusters == 2){
        #Need to sum up entries in the vector for number of clusters
        matrixnrow <- 0
        for (i in 1:input$waves){
          eval(parse(text=paste0("matrixnrow <- matrixnrow + input$waves.r1c",i)))
        }
      }
      #build empty matrix
      flexn <- matrix(NA,ncol=matrixncol,nrow=matrixnrow)
      #fill it in with HTML cells
      for (i in 1:matrixncol){#For each column i...
        for (j in 1:matrixnrow){#For each row j...
          eval(parse(text=paste0("flexn[", j, ",",i,"] <- input$Ntimecluster.r",j,"c",i)))
        }
      }
    }
    
    obj <- swPwr.FlexN(design=design,distn=input$distn,n=flexn,mu0=input$mu0,mu1=input$mu1,tau=input.tau,eta=input.eta,sigma=input$sigma,gamma=input.gamma,rho=input.rho,alpha=input$alpha)
    
    paste("Power:",round(obj,6))
  })
  
  ####################################################
  ##Other descriptive outputs
  ####################################################
  
  output$textoutput2 <-  renderText({
    vector.temp <- variance_components_vector()
    input.tau <- vector.temp[1]
    input.gamma <- vector.temp[2]
    input.eta <- vector.temp[3]
    input.rho <- vector.temp[4]
    
    if (input$variance_param == 1){
    if (input$eta == 0 & input$rho == 0){
      a <- paste("within-period ICC:",(input$tau^2+input$gamma^2)/(input$tau^2+input$gamma^2+input$sigma^2))
      b <- paste("CAC:",(input$tau^2)/(input$tau^2+input$gamma^2))
      
      if (input$icc == 0){
        b <- paste("CAC is not well-defined when ICC is 0; power depends only on ICC and sigma.")
      }
    }else{
      a <- paste("ICC and CAC are not well-defined for this design.")
      b <- paste("")
    }}
    
    if (input$variance_param == 2){
      a <- paste("within-period ICC:",input$icc)
      b <- paste("CAC:",input$cac)
    }
    
    f <- paste("Tau, Gamma, Eta, and Rho:",input.tau,input.gamma,input.eta,input.rho)
    
    #Total number of clusters
    if (input$clusters == 1){
      c <- paste("Total number of clusters:",input$nclust*input$waves)
    }
    #Total number of individuals
    if(input$whichN == 1 & input$clusters == 1){
      d <- paste("Total number of individuals (cross-sectional):",input$nclust*input$waves*(input$waves+input$extra.time+as.numeric(as.logical(input$all.ctl.time0)))*input$n)
    }
    if(input$clusters != 1 | input$whichN != 1){
      
      ##################Obnoxious matrix unpacking, copied from power calculation section and added definition of d within statements.
      
      tx.effect.num <- as.numeric(unlist(strsplit(input$tx.effect,",")))
      design.vector <- rep(1,length=input$waves)*input$nclust#This is the design vector for the simple case
      #Design vector for the case where the number of clusters varies by wave:
      if (input$clusters == 2){
        for (i in 1:input$waves){
          eval(parse(text=paste0("design.vector[i] <- input$waves.r1c",i)))
        }
      }
      
      design <- swDsn(clusters=design.vector,  extra.time=input$extra.time, all.ctl.time0=as.logical(input$all.ctl.time0), tx.effect=tx.effect.num)
      
      #Assign n as a scalar, vector, or matrix
      #Take input about the number of observations, and return a scalar, vector, or matrix.
      if (input$whichN == 1){#scalar
        flexn <- input$n
        d <- paste("Total number of individuals (cross-sectional):",sum(design.vector)*(input$waves+input$extra.time+as.numeric(as.logical(input$all.ctl.time0)))*input$n)
        
      } else if(input$whichN == 2){#vector (varies by cluster)
        #Identify dimension of vector
        #calculate number of clusters: depends on 'clusters'
        if (input$clusters == 1){
          vectorlength <- input$nclust*input$waves#number of clusters = number of clusters per wave * number of waves
        }
        if (input$clusters == 2){
          #Need to sum up entries in the vector for number of clusters
          vectorlength <- 0
          for (i in 1:input$waves){
            eval(parse(text=paste0("vectorlength <- vectorlength + input$waves.r1c",i)))
          }
        }
        #build empty vector
        flexn <- rep(NA,vectorlength)
        #fill it with HTML cells
        for (i in 1:vectorlength){
          eval(parse(text=paste0("flexn[", i, "] <- input$Ncluster.r",i,"c1")))
        }
        
        #Case where n varies by cluster, total number of individuals is:
        d <- paste("Total number of individuals (cross-sectional):",sum(flexn)*(input$waves+input$extra.time+as.numeric(as.logical(input$all.ctl.time0))))
      } else if(input$whichN == 3){#matrix
        #Identify dimension of matrix
        #calculate total time: number of waves + extra.time + all.ctl.time0
        matrixncol <- input$waves+input$extra.time+as.numeric(as.logical(input$all.ctl.time0))
        #calculate number of clusters: depends on 'clusters'
        if (input$clusters == 1){
          matrixnrow <- input$nclust*input$waves#number of clusters = number of clusters per wave * number of waves
        }
        if (input$clusters == 2){
          #Need to sum up entries in the vector for number of clusters
          matrixnrow <- 0
          for (i in 1:input$waves){
            eval(parse(text=paste0("matrixnrow <- matrixnrow + input$waves.r1c",i)))
          }
        }
        #build empty matrix
        flexn <- matrix(NA,ncol=matrixncol,nrow=matrixnrow)
        #fill it in with HTML cells
        for (i in 1:matrixncol){#For each column i...
          for (j in 1:matrixnrow){#For each row j...
            eval(parse(text=paste0("flexn[", j, ",",i,"] <- input$Ntimecluster.r",j,"c",i)))
          }
        }
        
        #Case where n varies by cluster and time, total number of individuals is:
        d <- paste("Total number of individuals (cross-sectional):",sum(flexn))
      }
      
      ###################
      if (input$clusters == 2){#number of clusters varies over time
        c <- paste("Total number of clusters:",sum(design.vector))
      }
      
    }
    
    e <- paste("Total number of time points:",input$waves+input$extra.time+as.numeric(as.logical(input$all.ctl.time0)))
    
    paste("Summary information (see help tab for details):",e,c,d,a,b,f,sep="\n")
    
  })
  
}

####################################################
##Create Shiny app
####################################################



shinyApp(ui = ui, server = server)
