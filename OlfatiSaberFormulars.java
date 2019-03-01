package de.dlr.swex.algorithm.saber;

import android.util.Log;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import de.dlr.swex.Vektor;
import de.dlr.swex.algorithm.IAgent;
import de.dlr.swex.gui.Obstacle;

/**
 * Olfati-Saber Flocking
 *
 * @see <a href="http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.121.7027&rep=rep1&type=pdf">http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.121.7027&rep=rep1&type=pdf</a>
 */

public class OlfatiSaberFormulars {


    public double SafetyBoxDistance = 1.0;
    public double MaxForce = 8.0f;
    public double MaxSpeed = 1.0f;

    public int a = 5;
    public float b = 5;
    public double c = Math.abs(a - b) / Math.sqrt(4 * a * b);
    public double h = 0.2f;
    public double EPSILON = 0.1f;
    public double d = 6;
    public double d_alpha = sigma_norm(d);
    public double r = 1.2 * d;
    public double r_alpha = sigma_norm(r);
    public float c1 = 5;
    public float c2 = 5;

  /*
       // Stable Parameters for not moving leader
    public int a = 5;
    public float b = 5;
    public double c = Math.abs(a - b) / Math.sqrt(4 * a * b);
    public double h = 0.5f;
    public double EPSILON = 0.9f;
    public double d = 4;
    public double d_alpha = sigma_norm(d);
    public double r = 1.05 * d;
    public double r_alpha = sigma_norm(r);
    public float c1 = 0.1f;
    public float c2 = 0; */


    // Addition Parameters for Algorithm3
    public float c1_alpha = 1;
    public float c1_beta = 5;
    public float c1_gamma = 5;
    public double c2_alpha = 2 * Math.sqrt(c1_alpha);
    public double c2_beta = 2 * Math.sqrt(c1_beta);
    public double c2_gamma = 1;

    public double d_straight = 2;
    public double d_beta = sigma_norm(d_straight);
    public double r_straight = 1.2 * d_straight;
    //public double r_beta = sigma_norm(r_straight);



    public OlfatiSaberFormulars() {

    }

    // Added for Algorithmus 3

    public double Phi_beta(double z) {
        return p_h(z / d_beta) * (sigma_1(z - d_beta) - 1);
    }

    public double b_ik(Vektor q) {
        return p_h(sigma_norm(q) / d_beta);
    }

    public Vektor sphericalObstaclePosition(Obstacle o, IAgent q_i) {
        Vektor v = Vektor.subtract(o.location, q_i.getPosition());
        double mu = o.radius / v.getLenght();
        Vektor position = q_i.getPosition();
        position.multiply(mu);

        Vektor second = new Vektor(o.location);
        second.multiply(1 - mu);

        position.add(second);
        return position;
    }

    public Vektor sphericalObstacleVelocity(Obstacle o, IAgent q_i) {
        Vektor p = q_i.getPosition();

        Vektor v = Vektor.subtract(o.location, p);
        double mu = 0;
        if (v.getLenght() != 0) {
            mu = o.radius / v.getLenght();
        }
        Vektor a_k = v;
        a_k.normalize();

        double[][] matrixData = {{a_k.getX()}, {a_k.getY()}, {a_k.getZ()}};
        RealMatrix ak = MatrixUtils.createRealMatrix(matrixData);
        RealMatrix akT = ak.transpose();
        RealMatrix identity = MatrixUtils.createRealIdentityMatrix(3);
        RealMatrix P = identity.subtract(ak.multiply(akT));

        Vektor velocity = new Vektor(q_i.getVelocity());
        velocity.multiply(mu);
        double[][] matrixData2 = {{velocity.getX()}, {velocity.getY()}, {velocity.getZ()}};
        RealMatrix velMatrix = MatrixUtils.createRealMatrix(matrixData2);
        RealMatrix result = P.multiply(velMatrix);

        RealVector vector = result.getColumnVector(0);
        double[] debug = vector.toArray();

        Vektor force = new Vektor(debug[0], debug[1], debug[2]);
        if (Double.isNaN(force.getX()) || Double.isNaN(force.getY()) || Double.isNaN(force.getZ())) {
            Log.wtf("dlr", "NaN Error 1");
            return new Vektor(0, 0, 0);
        }
        if (force.getX() != force.getX()) {
            Log.wtf("dlr", "NaN Error 2");
            return new Vektor(0, 0, 0);
        }
        return force;

    }
    // End Algorithm3 Formulars

    public double Phi_alpha(double z) {
        return p_h(z / r_alpha) * Phi(z - d_alpha);
    }


    private double p_h(double z) {
        if (z >= 0 && z < h) {
            return 1;
        } else if (z >= h && z <= 1) {
            return 0.5 * (1 + Math.cos(Math.PI * ((z - h) / (1 - h))));
        } else {
            return 0;
        }
    }


    private double Phi(double z) {
        return 0.5 * ((a + b) * sigma_1(z + c) + (a - b));
    }

    private double sigma_1(double z) {
        return z / Math.sqrt(1 + z * z);
    }

    public Vektor sigma_1ForAlgorithm3(Vektor z) {
        z.multiply(1 / Math.sqrt(1 + z.getLenght() * z.getLenght()));
        return z;
    }

    public double sigma_norm(Vektor v) {
        return (1 / EPSILON) * (Math.sqrt(1 + EPSILON * v.getLenght() * v.getLenght()) - 1);
    }

    public double sigma_norm(double z) {
        return sigma_norm(new Vektor(z, 0, 0));
    }

    public Vektor n_ij(Vektor v) {
        v.multiply(1 / Math.sqrt(1 + EPSILON * v.getLenght() * v.getLenght()));
        return v;
    }

    public double a_ij(Vektor agentToNeighbour) {
        return p_h(sigma_norm(agentToNeighbour) / r_alpha);
    }


}
