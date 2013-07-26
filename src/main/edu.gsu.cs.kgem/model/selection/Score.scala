package edu.gsu.cs.kgem.model.selection

import edu.gsu.cs.kgem.model.KGEM

/**
 * Created with IntelliJ IDEA.
 * User: nicholas
 * Date: 7/17/13
 * Time: 10:19 AM
 * To change this template use File | Settings | File Templates.
 */
object Score {
  /**
   * Scoring function type to make the model selection
   * higher-order parameter cleaner. After instantiation
   * it can be treated like a function
   * e.g. myScoringFunc(k, n)
   */

  type ScoringFunction =  (Int, Int) => Double

  /**
  * Constants for parsing.
  */
  val AIC_     = "AIC"
  val AICc_    = "AICc"
  val BIC_     = "BIC"
  val BICMAP_  = "BICMAP"

  /**
   * Get the scoring function by name.
   * Implemented functions are:
   *  AIC, AICc, BIC, BICMAP.
   * @param name
   *  The name of the scoring function to return.
   * @return
   *  The scoring function.
   */
  def getScoringFunction(name: String): ScoringFunction = name match {
    case AIC_ => AIC
    case AICc_ => AICc
    case BIC_ => BIC
    case BICMAP_ => BICMAP
    case _ => throw new IllegalArgumentException("Unknown scoring function!")
  }

  def getAlpha(name: String, isConservative: Boolean): Double = name match {
    case AIC_ | AICc_ | BIC_  => 0d
    case BICMAP_ => if (isConservative) 0.5 else 0.1
    case _ => throw new IllegalArgumentException("Unknown scoring function!")
  }
  /**
   * Compute the Akaike Information Criterion score:
   * 2(k - ln L)
   * @param k
   * Number of parameters to the model
   * @param n
   * Number of data-points.
   * @return
   * AIC score
   */
  def AIC(k: Int, n: Int): Double = {
    2.0 * (k - KGEM.loglikelihood)
  }

  /**
   * Compute the corrected Akaike Information Criterion score.
   * This penalizes harder than AIC when either n is small or
   * k is large.
   * 2(k - ln L) + (2k(k+1))/(n - k + 1)
   * @param k
   * Number of parameters to the model
   * @param n
   * Number of data-points
   * @return
   * Corrected AIC score
   */
  def AICc(k: Int, n: Int): Double = {
    val kk = k.toDouble
    val nn = n.toDouble
    2.0 * (k - KGEM.loglikelihood) + ((2 * kk * (kk + 1)) / (nn - kk + 1))
  }


  /**
   * Compute the Bayesian Information Criterion score.
   * This penalizes harder than AIC.
   * k ln n - 2 ln L
   * @param k
   * Number of parameters to the model
   * @param n
   * Number of data-points
   * @return
   * BIC score
   */
  def BIC(k: Int, n: Int): Double = {
    -2.0 * KGEM.loglikelihood + (k * math.log(n.toDouble))
  }

  /**
   * Compute the Bayesian Information Criterion score with
   * the MAP estimate instead of log Likelihood. Like BIC,
   * this penalizes harder than AIC.
   * k ln n - 2 ln MAP
   * @param k
   * Number of parameters to the model
   * @param n
   * Number of data-points
   * @return
   * BIC score used with MAP estimate
   */
  def BICMAP(k: Int, n: Int): Double = {
    -2.0 * KGEM.maximumAP + (k * math.log(n.toDouble))
  }
}
