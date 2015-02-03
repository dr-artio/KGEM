package edu.gsu.cs.kgem.model

/**
 * Created by Alex on 12/1/2014.
 */
trait SeedFinder {
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
  protected def hammingDistance(r1: Read, r2: Read): Double = {
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
  def hammingDistance(s: String, t: String): Double = {
    val l = s.length
    if (l != t.length) {
      throw new IllegalArgumentException("Hamming Distance: Strings have different lengths")
    }
    var r = 0.0
    for (i <- 0 until l) {
      if (s(i) != t(i) && s(i) != 'N' && t(i) != 'N' && t(i) != ' ' && s(i) != ' ') {
        r += 1
      } else if (s(i) != t(i) && (s(i) == 'N' || t(i) == 'N' || t(i) == ' ' || s(i) == ' '))  {
        //r += 0.2
      }
    }
    r
  }
}
