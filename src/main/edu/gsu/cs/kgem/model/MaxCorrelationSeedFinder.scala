package edu.gsu.cs.kgem.model

import java.io.File

import edu.gsu.cs.kgem.exec._
import org.biojava3.core.sequence.DNASequence
import org.biojava3.core.sequence.io.FastaWriterHelper

import scala.collection.mutable.ListBuffer
import scala.collection.JavaConverters._

/**
 * Author: Alexander Atryomenko <aartyomenko@cs.gsu.edu>
 * Initialization of candidates with the use of
 * Linkage Disequilibrium. Finds chain of correlated SNPs
 * to identify haplotypes.
 */
object MaxCorrelationSeedFinder extends SeedFinder {
  def findSeeds(reads: Iterable[Read], k: Int, threshold: Int, count_threshold: Int): Iterable[Genotype] = {
    val seeds = new ListBuffer[Genotype]()

    log("Initialize genotype...")
    val first = new Genotype(reads)
    log("Initial genotype constructed.")
    log("%s_eps: %f_size: %f".format(first.ID, first.epsilon, first.size))
    var i = 0
//    for (d <- first.data) {
//      val m = d.values.max
//      println("%d\t%.0f".format(i, d.values.sum - m))
//      i += 1
//    }

    seeds += first
    var cond = true

    while (seeds.size < k && !seeds.forall(_.convergen)) {
      log("Current size: %d".format(seeds.size))
      cond = false
      for (t <- seeds.filterNot(_.convergen)) {

        var newCenter = t.getSecondGenotype
        var secondGenotype: Genotype = null

        while (newCenter != None) {
          log("Step of second haplotype search")
          secondGenotype = newCenter.asInstanceOf[Genotype]
          newCenter = secondGenotype.getSecondGenotype
        }
        t.convergen = secondGenotype == null
        if (secondGenotype != null ) {
          for (r <- secondGenotype.reads) t.removeRead(r)
          val oldCenter = t.toIntegralString
          val newCenter = secondGenotype.toIntegralString
          log("SNPs: %f".format(hammingDistance(oldCenter, newCenter)))

          val newSeed = t.reads.par.filter(r => hammingDistance(newCenter, r.seq) < hammingDistance(oldCenter, r.seq)).seq
          log("New seed size: %f".format(newSeed.map(_.freq).sum))
          log("Initialize new seed...")
          if (newSeed.size > 1) {
            val seed = new Genotype(newSeed)
            log("initialized.")
            log("Update previous seed...")
            for (r <- newSeed)
              t.removeRead(r)
            log("Updated")

            log("New voronoi region size: %f".format(seed.size))
            seeds += seed
//            for (s <- seeds) {
//              log("id: %s eps: %f size: %f".format(s.ID, s.epsilon, s.size))
//              val records = s.reads.map(x => {
//                x.ids.map(t => {
//                  val rec = new DNASequence(x.seq)
//                  rec.setOriginalHeader(t)
//                  rec
//                })
//              }).flatten.asJava
//              //FastaWriterHelper.writeNucleotideSequence(new File("/research_data/sasha/IAV_UCLA/%s.fas".format(s.ID)), records)
//            }
          }
        }

      }


      Genotype.eps = seeds.map(g => g.epsilon * g.size).sum / seeds.map(_.size).sum

    }
    seeds.toList
  }
}