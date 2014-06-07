package edu.gsu.cs.kgem.model

import collection.mutable
import util.Random
import org.apache.commons.math3.distribution.BinomialDistribution
import edu.gsu.cs.kgem.exec.log

/**
 * Created with IntelliJ IDEA.
 * User: aartyomenko
 * Date: 3/23/13
 * Time: 8:26 PM
 * To change this template use File | Settings | File Templates.
 */
object KGEM {
  private var reads: List[Read] = null
  private var zreads: List[(Read, Int)] = null
  private var em: EM = null
  private var tr = 0.0
  private val pValue = 0.05
  var table: mutable.MutableList[Map[String, Iterable[(Read, Int)]]] = null
  var threshold = 0

  private def initSeeds(n: Int): List[Genotype] = {
    val seeds = sample(reads, n)
    seeds.map(s => new Genotype(s.seq)).toList
  }

  def initReads(reads: List[Read]) = {
    this.reads = reads
    zreads = reads.zipWithIndex
    this.table = new mutable.MutableList[Map[String, Iterable[(Read, Int)]]]()
    val ml = reads.map(r => r.end).max
    var i = 0
    val nucls = Genotype.sMap.keys
    while (i < ml) {
      table += nucls.map(nucl =>
        (nucl, zreads.filter(r => r._1.seq(i).equals(nucl(0))))).toMap
      i += 1
    }
  }

  def run(gens: Iterable[Genotype]) = {
    var collapse = gens
    var collapsed = gens.size
    var genMap = new mutable.HashMap[Genotype, Int]()
    do {
      collapsed = collapse.size
      runKgem(collapse)
      val col = collapse.map(g => (g.toIntegralString, g)).toMap
      collapse = col.values.toList
      collapse = thresholdClean(collapse, tr)
      collapse.view.zipWithIndex.foreach(pair => genMap += pair)
      collapsed -= collapse.size
      log("KGEM collapsed %d genotypes".format(collapsed))
    } while (collapsed > 0)
    collapse
  }

  def initThreshold(tr: Double) {
    this.tr = tr
    this.threshold = (reads.map(r => r.freq).sum * tr).toInt
    log("Set threshold: %f".format(tr))
  }

  def initThreshold(length: Int) = {
    this.tr = getThreshold(length)
    log("Computed threshold: %f".format(tr))
  }

  private def thresholdClean(gens: Iterable[Genotype], tr: Double) = {
    val cleaned = gens.filter(g => g.freq >= tr)
    cleaned
  }

  private def runKgem(gens: Iterable[Genotype]) = {
    for (g <- gens) g.convergen = false
    var i = 1
    while (!gens.forall(g => g.convergen) && i <= 5) {
      val st = System.currentTimeMillis
      rounding(gens)
      runEM(gens)
      alleleFreqEstimation(gens)
      log("KGEM iteration #%d done in %.2f minutes".format(i, (System.currentTimeMillis - st) * 1.0 / 60000))
      i += 1
    }
    rounding(gens)
  }

  private def rounding(gens: Iterable[Genotype]) = {
    for (g <- gens.par) g.round
  }

  private def runEM(gens: Iterable[Genotype]) = {
    em = new EM(gens.toList, reads.toList)
    em.run
  }

  private def alleleFreqEstimation(gens: Iterable[Genotype]) = {
    val pqrs = em.eStep
    val pargens = gens.view.zipWithIndex.par
    for (g <- pargens) doAlleleFreqEstimation(g._1, pqrs(g._2))
  }

  private def doAlleleFreqEstimation(g: Genotype, pqs: Array[Double]): Unit = {
    if (g.convergen) return
    val prev = g.toIntegralString
    for (e <- g.data.zipWithIndex) {
      for (v <- e._1.keys) e._1(v) = 0
      val idx = e._2
      val rs = table(idx)
      for (n <- rs.keys) {
        val res = rs(n)
        for (r <- res)
          e._1(n) += pqs(r._2)
      }
    }
    g.convergen = prev.equals(g.toIntegralString) || (g.freq < tr)
  }

  /**
   * Finds argmin_x((1-p_i)**L > 1 - pValue), where
   * p_i = F(x) is survival function of binomial distribution
   * pValue - measure of randomness (Default: 5% hardcoded)
   * n - number of reads in a given sample
   * @return
   * Estimated threshold
   */
  private def getThreshold(L: Int) = {
    val n = reads.view.map(_.freq).sum
    val topBound = 1 - pValue
    val p = Genotype.eps
    var step = (n / 2).toInt
    var x = step
    while (step > 1) {
      step /= 2
      val p_i = sf(x, n.toInt, p)
      if (Math.pow(1 - p_i, L) > topBound) x -= step
      else x += step
    }
    threshold = x
    x / n
  }

  /**
   * Method computing Survival Function
   * (1 - CDF(x)) for binomial distribution
   * with a given parameters
   * @param x
   * Argument for Pr(X>=x)
   * @param n
   * Number of samples in Binomial
   * distribution
   * @param p
   * Probability of success
   */
  private def sf(x: Int, n: Int, p: Double) = {
    1.0 - new BinomialDistribution(n, p).cumulativeProbability(x)
  }

  /**
   * Method for choosing random sample of size @size from the
   * collection of objects. Returns the whole list is @size
   * is greater than size of the original collection.
   * @param iter
   * Any iterable collection
   * @param size
   * Size of the required sample
   * @tparam T
   * Generic parameter (Class of objects in list)
   * @return
   * List of randomly chosen elements from collection
   * of specified size @size
   */
  private def sample[T](iter: Iterable[T], size: Int) = {
    if (iter.size < size) iter.toList
    var res = new mutable.MutableList[T]()
    val rnd = new Random(System.currentTimeMillis)
    var needed = size
    var len = iter.size
    val iterator = iter.iterator
    while (needed > 0 && iterator.hasNext) {
      val item = iterator.next()
      if (rnd.nextInt(len) < needed) {
        res += item
        needed -= 1
      }
      len -= 1
    }
    res
  }
}
