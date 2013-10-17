package edu.gsu.cs.kgem.model

import collection.mutable
import edu.gsu.cs.kgem.model.initialization.MaxDistanceSeedFinder
import util.Random
import org.apache.commons.math3.distribution.BinomialDistribution
import edu.gsu.cs.kgem.model.estimation.{EMMAP, EM}
import scala.collection.mutable.ListBuffer

/**
 * Created with IntelliJ IDEA.
 * User: aartyomenko
 * Date: 3/23/13
 * Time: 8:26 PM
 * To change this template use File | Settings | File Templates.
 */
object KGEM {
  private var reads: List[Read] = List[Read]()
  private var zreads = reads.zipWithIndex
  private var em = new EM(List[Genotype](), List[Read]())
  private var tr = 0.0005
  private val pValue = 0.05
  var table = new mutable.MutableList[Map[String, Iterable[(Read, Int)]]]()
  var loglikelihood = 0.0
  var maximumAP = 0.0

  def initSeeds(n: Int): List[Genotype] = {
    val seeds = sample(reads, n)
    return seeds.map(s => new Genotype(s.seq)).toList
  }

  def initReads(reads: List[Read]) = {
    this.reads = reads
    zreads = reads.zipWithIndex
    val ml = reads.map(r => r.end).max
    var i = 0
    val nucls = Genotype.sMap.keys
    while (i < ml) {
      table += ((for (nucl <- nucls)
      yield (nucl, zreads.filter(r => r._1.seq(i).equals(nucl(0))))).toMap)
      i += 1
    }
  }

  def run(gens: Iterable[Genotype], alpha: Double = 0) = {
    var collapse = gens
    var collapsed = gens.size
    var genMap = new mutable.HashMap[Genotype, Int]()
    do {
      collapsed = collapse.size
      runKgem(collapse, alpha)
      val col = collapse.map(g => (g.toIntegralString, g)).toMap
      collapse = col.values.toList
      collapse = thresholdClean(collapse, tr)
      collapse.zipWithIndex.foreach(pair => genMap += pair)
      collapsed -= collapse.size
      println("KGEM collapsed %d genotypes".format(collapsed))
    } while (collapsed > 0)
    calcLogLikelihood(collapse, genMap, alpha)
    collapse
  }

  def runCl(gens: Iterable[Genotype], k: Int, alpha: Double = 0) = {
    var clusters = run(gens, alpha)

    do {
      val bg = getBadGenotype(clusters)
      clusters = clusters.filter(c => c != bg)
      clusters = run(clusters, alpha)
    } while (clusters.size > k)
    clusters
  }

  def getBadGenotype(gens: Iterable[Genotype]) = {
    val pairs = new ListBuffer[(Genotype, Genotype)]()
    var gg = gens
    println("------")
    while(!gg.tail.isEmpty){
      pairs ++= gg.tail.map(g => (gg.head, g))
      gg = gg.tail
      println(gg.tail.size)
      println(gg.head.ID)
    }
    println(pairs.size)
    val pair = pairs.toList.minBy(p => {
      val d = MaxDistanceSeedFinder.hammingDistance(p._1.toIntegralString, p._2.toIntegralString)*Math.sqrt(p._1.freq * p._2.freq)
      println("%.5f %d %d".format(d, p._1.ID, p._2.ID))
      d
    })
    if (pair._1.freq > pair._2.freq) pair._2
    else pair._1
  }

  def initThreshold(tr: Double) {
    this.tr = tr
    println("Set threshold: %f".format(tr))
  }

  def initThreshold = {
    this.tr = getThreshold
    println("Computed threshold: %f".format(tr))
  }

  private def calcLogLikelihood(gens: Iterable[Genotype], genIdxMap: mutable.HashMap[Genotype, Int], alpha: Double) = {
    val pqrs = em.eStep
    var ll = 0.0
    for (readpair <- zreads) {
      val read = readpair._1
      val idx = readpair._2
      ll += (read.freq * math.log(gens.map(gen => gen.freq * pqrs(genIdxMap(gen))(idx)).sum))
    }
    if (alpha == 0.0) {
      this.loglikelihood = ll
    } else {
      this.maximumAP = ll
    }
  }

  private def thresholdClean(gens: Iterable[Genotype], tr: Double) = {
    val cleaned = gens.filter(g => g.freq >= tr)
    cleaned
  }

  private def runKgem(gens: Iterable[Genotype], alpha: Double) = {
    for (g <- gens) g.convergen = false
    var i = 1
    while (!gens.forall(g => g.convergen) && i <= 10) {
      val st = System.currentTimeMillis
      rounding(gens)
      if (alpha > 0) runEM(gens, alpha)
        else runEM(gens)
      alleleFreqEstimation(gens)
      println("KGEM iteration #%d done in %.2f minutes".format(i, ((System.currentTimeMillis - st) * 1.0 / 60000)))
      i += 1
    }
  }

  private def rounding(gens: Iterable[Genotype]) = {
    for (g <- gens) g.round
  }

  def runEM(gens: Iterable[Genotype], alpha: Double) = {
    em = new EMMAP(gens.toList, reads.toList, alpha)
    em.run
  }

  def runEM(gens: Iterable[Genotype]) = {
    em = new EM(gens.toList, reads)
    em.run
  }

  private def alleleFreqEstimation(gens: Iterable[Genotype]) = {
    val pqrs = em.eStep
    for (g <- gens.zipWithIndex.par) doAlleleFreqEstimation(g._1, pqrs(g._2))
  }

  private def doAlleleFreqEstimation(g: Genotype, pqs: Array[Double]): Unit = {
    if (g.convergen) return
    val prev = g.toIntegralString
    for (e <- g.data.zipWithIndex.par) {
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
   * Finds argmin_x((F(x) < pValue/n), where
   * F(x) is survival function of binomial distribution
   * pValue - measure of randomness (Default: 5% hardcoded)
   * n - number of reads in a given sample
   * @return
   *         Estimated threshold
   */
  private def getThreshold = {
    val n = reads.map(r => r.freq).sum
    val topBound = pValue / n
    val p = em.eps
    var step = (n/2).toInt
    var x = step
    while (step > 1) {
      step /= 2
      if (sf(x, n.toInt, p)<topBound) x -= step
      else x += step
    }
    x / n
  }

  /**
   * Method computing Survival Function
   * (1 - CDF(x)) for binomial distribution
   * with a given parameters
   * @param x
   *          Argument for Pr(X>=x)
   * @param n
   *          Number of samples in Binomial
   *          distribution
   * @param p
   *          Probability of success
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
      val item = iterator.next
      if (rnd.nextInt(len) < needed) {
        res += item
        needed -= 1
      }
      len -= 1
    }
    res
  }

  /**
   * Get pqrs for generating clustering string
   * @return
   *         Matrix with P_qr 's
   */
  def getPqrs = {
    em.eStep
  }

  def getReads = {
    reads
  }
}
