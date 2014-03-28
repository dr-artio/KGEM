package edu.gsu.cs.kgem.model

import MaxDistanceSeedFinder.hammingDistance
import com.apporiented.algorithm.clustering.{Cluster, DefaultClusteringAlgorithm}
import scala.collection.mutable.ListBuffer
import org.biojava3.core.sequence.DNASequence
import org.apache.commons.lang3.StringUtils._
import edu.gsu.cs.kgem.exec.log
import collection.JavaConversions._

/**
 * Created by aartsiomenka1 on 3/28/2014.
 */
object Clusterer {
  val AMP = "&"
  /**
   * Perform Ward's clustering on obtained
   * cluster centers to fold cluster set to
   * desired size
   * @param genSeqs
   * Current cluster centers
   * @param sizes
   * Sizes of original clusters
   * @param k
   * Desired size of clusters
   * @return
   * Map with ids of original cluster
   * centers' ids and new id as a value
   */
  private def hierarchicalClustering(genSeqs: List[(Genotype, String)], sizes: Array[Int], k: Int) = {
    val l = genSeqs.length
    val distanceMatrix = Array.tabulate[Double](l, l)((i, j) => hammingDistance(genSeqs(i)._1.toIntegralString, genSeqs(j)._1.toIntegralString))

    val alg = new DefaultClusteringAlgorithm()

    val names = genSeqs.map(_._1.ID.toString).toArray[String]
    val cluster = alg.performClustering(distanceMatrix, sizes, names)

    var clusters = new ListBuffer[Cluster]()
    clusters += cluster
    while (clusters.size < k) {
      val nextSplitCluster = clusters.filter(_.getDistance != null).maxBy(_.getDistance)
      clusters.-=(nextSplitCluster)
      clusters ++= nextSplitCluster.getChildren
    }

    genSeqs.map(g => (g._1.ID, clusters.indexWhere(_.getName.split(AMP).exists(_.equals(g._1.ID.toString))))).toMap
  }

  /**
   * Cluster original reads
   * @param gens
   * Collection of Haplotypes
   * @param reads
   * Collection of reads
   * @param pqrs
   * Clustering coefficients
   */
  def getClusteredFastaSequencies(gens: Iterable[Genotype], reads: Iterable[Read],
                           pqrs: Array[Array[Double]], rs: Map[String, DNASequence], k: Int) = {
    val genSeqs = gens.map(g => (g, g.toIntegralString.replaceAll("-", "").replaceAll("N", ""))).toList

    val groups = cluster(reads, pqrs, genSeqs, rs)

    val sizes = getSizesForWardsClustering(groups, genSeqs)

    val clusterMap = hierarchicalClustering(genSeqs, sizes, k)

    val sortedNgReads = groups.groupBy(x => clusterMap(x._1.ID))

    sortedNgReads.map(cd => {
      cd._2.map(rd => {
        rd._2.map(r => {
          val ds = rs(r._1.ids.head)
          ds.setOriginalHeader("c%d_h%d_%s_%d".format(clusterMap(rd._1.ID), rd._1.ID, ds.getOriginalHeader, r._1.ids.size))
          ds
        })
      })
    }).flatten.flatten
  }

  private def getSizesForWardsClustering(groups: Map[Genotype, Iterable[(Read, Int)]],
                                     genSeqs: List[(Genotype, String)]) = {
    val sizes = new Array[Int](genSeqs.length)
    groups.foreach(g => sizes(genSeqs.indexWhere(t => t._1.ID == g._1.ID)) = g._2.map(_._1.freq).sum.toInt)
    log("Total number of reads: %d".format(sizes.sum))
    sizes
  }

  private def cluster(reads: Iterable[Read], pqrs: Array[Array[Double]],
                      genSeqs: List[(Genotype, String)],
                      rs: Map[String, DNASequence]) = {
    val groups = group(reads, pqrs, genSeqs)

    val dsThresholds = getThresholds(groups, genSeqs, rs)

    val bReads = getSuspeciousReads(groups, genSeqs, dsThresholds, rs)

    val gbReads = fixSuspeciousClusteredReads(bReads, genSeqs, dsThresholds, rs)

    val res =  finalizeClustering(reads, pqrs, gbReads, genSeqs, rs)

    res
  }

  private def getSuspeciousReads(groups: Map[Genotype, Iterable[(Read, Int)]],
                                 genSeqs: List[(Genotype, String)],
                                 dsThresholds: Map[Genotype, Double],
                                 rs: Map[String, DNASequence]) = {
    groups.map(g => {
      val threshold = dsThresholds(g._1)
      def readSeq(r: Read) = rs(r.ids.head).getSequenceAsString
      val genSeq = genSeqs.find(_._1 == g._1).head._2
      g._2.filter(r => getLevenshteinDistance(readSeq(r._1), genSeq) > threshold)
    }).flatten
  }

  private def getThresholds(groups: Map[Genotype, Iterable[(Read, Int)]],
                            genSeqs: List[(Genotype, String)],
                            rs: Map[String, DNASequence]) = {
    groups.map(g => {
      def readSeq(r: Read) = rs(r.ids.head).getSequenceAsString
      val genSeq = genSeqs.find(_._1 == g._1).head._2
      val dses = g._2.map(r => getLevenshteinDistance(readSeq(r._1), genSeq))
      val threshold = mean(dses) + 3 * sigma(dses)
      (g._1, threshold)
    }).toMap
  }

  private def finalizeClustering(reads: Iterable[Read], pqrs: Array[Array[Double]],
                                 gbReads: Map[(Read, Int), Genotype],
                                 genSeqs: List[(Genotype, String)],
                                 rs: Map[String, DNASequence]) = {
    reads.zipWithIndex.filter(x => !gbReads.contains(x) || gbReads(x) != null).groupBy(rd => {
      if (gbReads.contains(rd)) {
        gbReads(rd)
      } else {
        val cl = column(pqrs, rd._2)
        val index = cl.indexWhere(p => p == cl.max)
        genSeqs(index)._1
      }
    }).toMap
  }

  private def fixSuspeciousClusteredReads(bReads: Iterable[(Read, Int)],
                                          genSeqs: List[(Genotype, String)],
                                          dsThresholds: Map[Genotype, Double],
                                          rs: Map[String, DNASequence]) = {
    log("Suspicious clustered reads:")
    bReads.map(r => {
      log(r._1.ids.head)
      def readSeq(r: Read) = rs(r.ids.head).getSequenceAsString
      val genOp = genSeqs.find(g => getLevenshteinDistance(g._2, readSeq(r._1)) < dsThresholds(g._1))
      if (genOp.size == 1)
        (r, genOp.head._1)
      else
        (r, null)
    }).toMap
  }

  private def group(reads: Iterable[Read], pqrs: Array[Array[Double]],
                    genSeqs: List[(Genotype, String)]): Map[Genotype, Iterable[(Read, Int)]] = {
    reads.zipWithIndex.groupBy(rd => {
      val cl = column(pqrs, rd._2)
      val index = cl.indexWhere(p => p == cl.max)
      genSeqs(index)._1
    }).toMap
  }

  private def column[A, M[_]](matrix: M[M[A]], colIdx: Int)
                             (implicit v1: M[M[A]] => Seq[M[A]], v2: M[A] => Seq[A]): Seq[A] =
    matrix.map(_(colIdx))


  /**
   * Merge two maps with collections as
   * values. Result map will contain union
   * of keys as new keyset and joint lists of
   * values
   * @param m1
   * Operand 1
   * @param m2
   * Operand 2
   * @tparam K
   * Generic parameter type of Key
   * @tparam V
   * Generic parameter type of Value
   * @return
   */
  def merge[K, V](m1: Map[K, Iterable[V]], m2: Map[K, Iterable[V]]): Map[K, Iterable[V]] = {
    val k1 = Set(m1.keysIterator.toList: _*)
    val k2 = Set(m2.keysIterator.toList: _*)
    val intersection = k1 & k2

    val r1 = for (key <- intersection) yield (key -> (m1(key) ++ m2(key)))
    val r2 = m1.filterKeys(!intersection.contains(_)) ++ m2.filterKeys(!intersection.contains(_))
    r2 ++ r1
  }

  def mean[T <% Double](vals: Iterable[T]) = {
    val N = vals.size
    vals.map(_.toDouble).sum / N
  }

  def sigma[T <% Double](vals: Iterable[T]) = {
    val meanValue = mean(vals)
    Math.sqrt(mean(vals.map(x => Math.pow(x - meanValue, 2))))
  }
}
