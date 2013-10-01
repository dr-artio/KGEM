package cc.mallet.util;
/* Copyright (C) 2002 Univ. of Massachusetts Amherst, Computer Science Dept.
   This file is part of "MALLET" (MAchine Learning for LanguagE Toolkit).
   http://www.cs.umass.edu/~mccallum/mallet
   This software is provided under the terms of the Common Public License,
   version 1.0, as published by http://www.opensource.org.  For further
   information, see the file `LICENSE' included with this distribution. */


/**
 @author Andrew McCallum <a href="mailto:mccallum@cs.umass.edu">mccallum@cs.umass.edu</a>
 */

public class Dirichlet extends java.util.Random {
    double magnitude = 1;
    double[] partition;

    public Dirichlet(int seed) {
        super(seed);
    }

    public Dirichlet() {
        super();
    }

    public Dirichlet(int size, double alpha)
    {
        super();
        magnitude = size * alpha;
        partition = new double[size];
        partition[0] = 1.0 / size;
        for (int i=1; i<size; i++) {
            partition[i] = partition[0];
        }
    }

    public double[] nextDistribution() {
        double distribution[] = new double[partition.length];
        //For each dimension, draw a sample from Gamma(mp_i, 1)
        double sum = 0;
        for (int i=0; i<distribution.length; i++) {
            distribution[i] = this.nextGamma(partition[i] * magnitude, 1);
            if (distribution[i] <= 0) {
                distribution[i] = 0.0001;
            }
            sum += distribution[i];
        }

        //Normalize
        for (int i=0; i<distribution.length; i++) {
            distribution[i] /= sum;
        }

        return distribution;
    }

    /** Return a random double in the range 0 to 1, inclusive, uniformly sampled from that range.
     * The mean of this distribution is 0.5.  The variance is 1/12. */
    public synchronized double nextUniform() {
        long l = ((long)(next(26)) << 27) + next(27);
        return l / (double)(1L << 53);
    }

    /** Return a random double in the range a to b, inclusive, uniformly sampled from that range.
     * The mean of this distribution is (b-a)/2.  The variance is (b-a)^2/12 */
    public synchronized double nextUniform(double a,double b) {
        return a + (b-a)*nextUniform();
    }

    // generate Gamma(1,1)
    // E(X)=1 ; Var(X)=1
    /** Return a random double drawn from a Gamma distribution with mean 1.0 and variance 1.0. */
    public synchronized double nextGamma() {
        return nextGamma(1,1,0);
    }

    /** Return a random double drawn from a Gamma distribution with mean alpha and variance 1.0. */
    public synchronized double nextGamma(double alpha) {
        return nextGamma(alpha,1,0);
    }

    /** Return a random double drawn from a Gamma distribution with mean alpha*beta and variance alpha*beta^2. */
    public synchronized double nextGamma(double alpha, double beta) {
        return nextGamma(alpha,beta,0);
    }

    /** Return a random double drawn from a Gamma distribution
     *  with mean alpha*beta+lamba and variance alpha*beta^2.
     *  Note that this means the pdf is:
     *     <code>frac{ x^{alpha-1} exp(-x/beta) }{ beta^alpha Gamma(alpha) }</code>
     *  in other words, beta is a "scale" parameter. An alternative
     *  parameterization would use 1/beta, the "rate" parameter.
     */
    public synchronized double nextGamma(double alpha, double beta, double lambda) {
        double gamma=0;
        if (alpha <= 0 || beta <= 0) {
            throw new IllegalArgumentException ("alpha and beta must be strictly positive.");
        }
        if (alpha < 1) {
            double b,p;
            boolean flag = false;

            b = 1 + alpha * Math.exp(-1);

            while (!flag) {
                p = b * nextUniform();
                if (p > 1) {
                    gamma = -Math.log((b - p) / alpha);
                    if (nextUniform() <= Math.pow(gamma, alpha - 1)) {
                        flag = true;
                    }
                }
                else {
                    gamma = Math.pow(p, 1.0/alpha);
                    if (nextUniform() <= Math.exp(-gamma)) {
                        flag = true;
                    }
                }
            }
        }
        else if (alpha == 1) {
            // Gamma(1) is equivalent to Exponential(1). We can
            //  sample from an exponential by inverting the CDF:

            gamma = -Math.log (nextUniform ());

            // There is no known closed form for Gamma(alpha != 1)...
        }
        else {

            // This is Best's algorithm: see pg 410 of
            //  Luc Devroye's "non-uniform random variate generation"
            // This algorithm is constant time for alpha > 1.

            double b = alpha - 1;
            double c = 3 * alpha - 0.75;

            double u, v;
            double w, y, z;

            boolean accept = false;

            while (! accept) {
                u = nextUniform();
                v = nextUniform();

                w = u * (1 - u);
                y = Math.sqrt( c / w ) * (u - 0.5);
                gamma = b + y;

                if (gamma >= 0.0) {
                    z = 64 * w * w * w * v * v;  // ie: 64 * w^3 v^2

                    accept = z <= 1.0 - ((2 * y * y) / gamma);

                    if (! accept) {
                        accept = (Math.log(z) <=
                                2 * (b * Math.log(gamma / b) - y));
                    }
                }
            }

			/* // Old version, uses time linear in alpha
			   double y = -Math.log (nextUniform ());
			   while (nextUniform () > Math.pow (y * Math.exp (1 - y), alpha - 1))
			   y = -Math.log (nextUniform ());
			   gamma = alpha * y;
			*/
        }
        return beta*gamma+lambda;
    }
}

