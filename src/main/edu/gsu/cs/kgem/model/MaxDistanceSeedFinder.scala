package edu.gsu.cs.kgem.model

import java.lang.Math.min

import edu.gsu.cs.kgem.exec.log

import scala.collection.mutable

/**
 * Created with IntelliJ IDEA.
 * User: aartyomenko
 * Date: 4/13/13
 * Time: 1:43 PM
 * Class to wrap derandomized KGEM, i. e. initialization with maximizing distance between
 * seeds and detection of size of the population via distance threshold.
 */
object MaxDistanceSeedFinder extends SeedFinder {

  /**
   * Find seeds according to maximization Hamming distance between all pairs
   * of reads with pre-specified threshold. This is a 2-approximation to the
   * metric k-center problem.
   *
   * @param reads
   * Collection of reads
   * @param k
   * Maximum size of sample
   * @param threshold
   * Min hamming distance between seeds
   * @return
   */
  def findSeeds(reads: Iterable[Read], k: Int, threshold: Int, count_threshold: Int): Iterable[Genotype] = {
    val readArr = reads.toArray
    val first = getFirstSeed(readArr)
    var seeds = new mutable.MutableList[Genotype]()
    seeds += new Genotype(first.seq)
    var distanceMap = readArr.view.filter(r => !r.equals(first)).map(r => r -> hammingDistance(first.seq, r.seq) * r.freq).toMap
    var maxHD = Double.MaxValue
    var count = maxHD
    log("Count threshold: %d".format(count_threshold))
    while (seeds.size < k && (maxHD >= threshold && count > count_threshold)) {
      if (distanceMap.isEmpty) return seeds//.map(r => new Genotype(r.seq))
      val cur = distanceMap.filter(x => x._1.seq.count(_ != 'N')>= x._1.len * 0.85).maxBy(r => r._2 * r._1.freq)
      maxHD = distanceMap.view.map(r => r._2 / r._1.freq).max

      //if (maxHD < threshold) count = readArr.view.filter(r => hammingDistance(r, cur._1) < distanceMap(r) / r.freq).map(_.freq).sum
      val voronoi = readArr.view.filter(r => hammingDistance(r, cur._1) < distanceMap(r) / r.freq).toList
      val center = new Genotype(voronoi)
      count = voronoi.map(_.freq).sum
      val seed = seeds.find(s => hammingDistance(s.toIntegralString, center.toIntegralString) == 0)
      log(seed.toString)
      if (seed == None) seeds += center
      log("Read to center dist: %.2f".format(hammingDistance(cur._1.seq, center.toIntegralString)))
      distanceMap = distanceMap.map(e => e._1 -> min(e._2, hammingDistance(cur._1, e._1) * e._1.freq)).toMap
    }

    log("Estimated k: %d".format(seeds.size))
    log("Final max distance: %.0f count: %.0f".format(maxHD, count))
    seeds//.map(r => new Genotype(r.seq))
  }

  /**
   * Select first read randomly
   * @param readArr
   * Array with all reads
   * @return
   * one read
   */
  private def getFirstSeed(readArr: Array[Read]) = {
    val mc = readArr.map(r => r.seq.count(_ != 'N')).max
    val candidates = readArr.filter(r => r.seq.count(_ != 'N') == mc)
    candidates.maxBy(_.freq)
    //val s = candidates.size
    //val rnd = new Random()
    //candidates(rnd.nextInt(s))
  }
}
