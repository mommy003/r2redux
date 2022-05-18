
  #' cc_trf function
  #'
  #' This function transforms the predictive ability and 
  #' its standard error between the observed scale 
  #' and liability scale cc_trf(R2,se, K,P), where 
  #' R2 is the R squared value of the model for either 
  #' observed scale or liability scale, se is the 
  #' standard error for either observed scale or liability scale, 
  #' K and P is the prevalence of cases in the population.
  #' 
  #' @keywords Transformation of R2 between observed scale and liability scale
  #' @references
  #' Lee, S. H., Goddard, M. E., Wray, N. R., and Visscher, P. M. A better coefficient of determination for genetic profile analysis. Genetic epidemiology,(2012). 36(3): p. 214-224.
  #' @export 
  #' @examples 
  #' To get the transformed predictive ability 
  #' output=cc_trf(0.06, 0.002, 0.05, 0.05)
  #'
  #' output$R2l (transformed R2 to liability scale, assuming 0.06 in observed scale)
  #' 0.2679337
  #' 
  #' output$sel (transformed se to liability scale, assuming 0.002 is observed scale)
  #' 0.008931123
  #'
  #' output$R2O (transformed R2 to observed scale, assuming 0.06 in liability scale)
  #' 0.01343616
  #'
  #' output$seO (transformed to in observed scale, assuming 0.002 is liability scale)
  #' 0.000447872
  #'
  

cc_trf = function (R2,se, K,P) {

    thd = -qnorm(K,0,1)
    zv = dnorm(thd) #z (normal density)
    mv = zv/K #mean liability for case
    mv2 = -mv*K/(1-K) #mean liability for controls


    #R2 on the observed scale
    theta = mv*(P-K)/(1-K)*(mv*(P-K)/(1-K)-thd) #theta in equation (15)
    cv = K*(1-K)/zv^2*K*(1-K)/(P*(1-P)) #C in equation (15)
    R2l = R2*cv/(1+R2*theta*cv)
    R2O = R2/(cv-R2*theta*cv)


    #SE on the liability (From a Taylor series expansion)
    #var(h2l_r2) = [d(h2l_r2)/d(R2v)]^2*var(R2v) with d being calculus differentiation
    #sel = cv*(1-R2*theta)*se
    sel = (cv/(cv*R2*theta+1)^2)*se
    seO = (1/(cv*(theta*R2-1)^2))*se

    z=list(R2l=R2l,sel=sel,R2O=R2O,seO=seO)
    return(z)

  }

