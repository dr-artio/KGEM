package edu.gsu.cs.kgem.model

import java.util.Random
import scala.collection.mutable

/**
 * Created with IntelliJ IDEA.
 * User: aartyomenko
 * Date: 4/13/13
 * Time: 1:43 PM
 * Class to wrap derandomized KGEM, i. e. initialization with maximizing distance between
 * seeds and detection of size of the population via distance threshold.
 */
object MaxDistanceWrapper {

  /**
   * Run KGEM procedure for the reads
   * @param reads
   * Collection of reads
   * @param k
   * Maximum size of population
   * @param threshold
   * Min hamming distance between seeds
   * @return
   * Collection of genotypes (haplotypes)
   */
  def run(reads: Iterable[Read], k: Int, threshold: Int) = {
    val gens = initSeeds(reads, k, threshold, null).map(r => new Genotype(r.seq))
    println("Initialization gave %d seeds in the population(Threshold: %d, k: %d)".format(gens.size, threshold, k))
    KGEM.run(gens.toList)
  }

  /**
   * Run KGEM procedure for the reads and consensus as first seed
   * @param reads
   * Collection of reads
   * @param k
   * Maximum size of population
   * @param threshold
   * Min hamming distance between seeds
   * @param consensus
   * Consensus wrapped into Read object
   * @return
   * Collection of genotypes (haplotypes)
   */
  def run(reads: Iterable[Read], k: Int, threshold: Int, seeds: Iterable[Read]) = {
    val gens = initSeeds(reads, k, threshold, seeds).map(r => new Genotype(r.seq))
    println("Initialization gave %d seeds in the population(Threshold: %d, k: %d)".format(gens.size, threshold, k))
    KGEM.run(gens.toList)
  }

  /**
   * Initialize seeds according to maximization Hamming
   * distance between all pairs of reads with pre-
   * specified threshold
   *
   * @param reads
   * Collection of reads
   * @param k
   * Maximum size of sample
   * @param threshold
   * Min hamming distance between seeds
   * @return
   */
  private def initSeeds(reads: Iterable[Read], k: Int, threshold: Int, candidates: Iterable[Read]): mutable.MutableList[Read] = {
    if (candidates != null) {
      val seeds = new mutable.MutableList[Read]()
      for (c <- candidates) seeds += c
      return seeds
    }
    val readArr = reads.toArray
    val first = getFirstSeed(readArr)
    var seeds = new mutable.MutableList[Read]()
    seeds += first
    var distanceMap = readArr.filter(r => !r.equals(first)).map(r => (r, hammingDistance(first, r)))
    while (seeds.size < k) {
      if (distanceMap.isEmpty) return seeds
      val cur = distanceMap.maxBy(e => e._2)
      println("Current max HD: %d".format(cur._2))
      if (cur._2 < threshold) return seeds
      seeds += cur._1
      distanceMap = distanceMap.map(e => (e._1, min(e._2, hammingDistance(cur._1, e._1)))) //.filter(e => (e._2 >= threshold))
      //println("Current size of DistnceMap: %d".format(distanceMap.size))
      //distanceMap = readArr.filter(r => !seeds.contains(r)).map(r => (r, seeds.map(s => hammingDistance(s, r)).min)).toMap
    }
    return seeds
  }

  /**
   * Select first read randomly
   * @param readArr
   * Array with all reads
   * @return
   * one read
   */
  private def getFirstSeed(readArr: Array[Read]) = {
    val s = readArr.size
    val rnd = new Random()
    readArr(rnd.nextInt(s))
  }

  @inline
  def min(i1: Int, i2: Int): Int = {
    if (i1 < i2) return i1
    i2
  }

  /**
   * Wrapper for hamming distance between reads
   * @param r1
   * Read 1
   * @param r2
   * Read 2
   * @return
   * Hamming Distance between reads
   */
  @inline
  private def hammingDistance(r1: Read, r2: Read): Int = {
    if (r1.equals(r2)) return 0
    hammingDistance(r1.seq, r2.seq)
  }

  /**
   * Compute hamming distance between two strings
   * of the same length
   * @param s
   * String 1
   * @param t
   * String 2
   * @return
   * Hamming distance between s and t if
   * their length is the same and -1
   * otherwise
   */
  @inline
  def hammingDistance(s: String, t: String): Int = {
    val l = s.length
    if (l != t.length) {
      System.err.println("Hamming Distance: Strings have different length")
      return -1
    }
    var r = 0
    for (i <- 0 until l) {
      if (s(i) != t(i) && s(i) != ' ' && t(i) != ' ' && t(i) != '-' && s(i) != '-') {
        r += 1
      }
    }
    return r
  }
}
