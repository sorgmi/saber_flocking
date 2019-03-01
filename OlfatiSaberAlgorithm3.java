package de.dlr.swex.algorithm.saber;

import java.util.ArrayList;

import de.dlr.swex.Vektor;
import de.dlr.swex.algorithm.IAgent;
import de.dlr.swex.algorithm.parameter.OlfatiSaberParameterDialogAlgorithm3;
import de.dlr.swex.gui.Obstacle;
import de.dlr.swex.ros.FloatPublisher;
import de.dlr.swex.ros.RosLoggerPublisher;
import rosgraph_msgs.Log;

/**
 * Olfati-Saber Flocking - Algorithm 3 with Obstacle Avoidance
 *
 * @see <a href="http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.121.7027&rep=rep1&type=pdf">http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.121.7027&rep=rep1&type=pdf</a>
 */

public class OlfatiSaberAlgorithm3 {


    private ArrayList<IAgent> agents = new ArrayList<>();
    private IAgent leader;
    private OlfatiSaberFormulars formulars = new OlfatiSaberFormulars();
    private ArrayList<Obstacle> obstacles = new ArrayList<>();

    private float HIGH_VALUE = 1000;
    private float AGENT_RADIUS = 0.25f;
    private FloatPublisher floatPublisherminDist = new FloatPublisher("swarmTool/debug/minDist");
    private FloatPublisher floatPublisherminZ = new FloatPublisher("swarmTool/debug/minZ");
    private FloatPublisher floatPublisherAvgLeader = new FloatPublisher("swarmTool/debug/avgDistLeader");
    private FloatPublisher floatPublisherminDistObstacle = new FloatPublisher("swarmTool/debug/minDistObstacle");


    public OlfatiSaberAlgorithm3(ArrayList<IAgent> list, ArrayList<Obstacle> obstacles) {
        this.agents = list;
        this.obstacles = obstacles;
    }


    public void iterate() {
        leader = null;
        for (IAgent q : agents) {
            if (q.isSelected()) leader = q;
        }
        for (IAgent agent : agents) {
            Vektor acceleration = agent.getAccelaration();
            Vektor velocity = agent.getVelocity();

            ArrayList<IAgent> neighbours = getNeighboursFrom(agent);
            Vektor fig = fig(agent, neighbours);
            Vektor fid = fid(agent, neighbours);
            fig.add(fid);


            Vektor u_ia = fig;
            Vektor u_iB = u_iB(agent, neighbours);
            Vektor u_iY = u_iY(agent);

            acceleration.add(u_ia);
            acceleration.add(u_iB);
            acceleration.add(u_iY);
            acceleration.limit(formulars.MaxForce);
            velocity.add(acceleration);
            velocity.limit(formulars.MaxSpeed);
            Vektor pos = agent.getPosition();
            pos.add(velocity);
            safetyBox(agent,pos);
            acceleration.multiply(0);
        }
        analyzeMin();
        analyzeAvgLeader();
        analyzeMinDistObst();
        //analyzeDeviationEnergy();
    }

    private void safetyBox(IAgent agent, Vektor pos){
        if (!agent.isSelected() && leader != null) {
           boolean error = false;
            for (IAgent neighbour : agents) {
                double distance = Vektor.subtract(agent.getPosition(), neighbour.getPosition()).getLenght();
                if (distance < formulars.SafetyBoxDistance && neighbour != agent) {
                    error = true;
                    agent.sendWaypoint((float)agent.getPosition().getX(),(float)agent.getPosition().getY(),(float)agent.getPosition().getZ(),0f,2f,2f);
                    if(neighbour != leader){
                        neighbour.sendWaypoint((float) neighbour.getPosition().getX(),(float) neighbour.getPosition().getY(),(float) neighbour.getPosition().getZ(),0f,2f,2f);
                    }
                }
            }
            if(!error){
                agent.sendWaypoint((float) pos.getX(), (float) pos.getY(), (float) leader.getPosition().getZ(), 0, (float) formulars.MaxSpeed, (float)formulars.MaxSpeed);
            }
            else{
                RosLoggerPublisher.getInstance().log("Safety Box ERROR", Log.FATAL);
            }

        }
    }

    /**
     * Follow Leader Force
     *
     * @param agent current Agent
     * @return Vektor force
     */
    private Vektor u_iY(IAgent agent) {
        if(leader != null){
            Vektor f1 = formulars.sigma_1ForAlgorithm3(Vektor.subtract(leader.getPosition(), agent.getPosition()));
            f1.multiply(-formulars.c1_gamma);

            Vektor f2 = Vektor.subtract(leader.getVelocity(), agent.getVelocity());
            f2.multiply(-formulars.c2_gamma);

            f1.add(f2);
            return f1;
        }
        else {
            return new Vektor(0,0,0);
        }
    }

    /**
     * Avoid Obstacles Force
     *
     * @param agent      current Agent
     * @param neighbours List with all Neighbours
     * @return Vektor force
     */
    private Vektor u_iB(IAgent agent, ArrayList<IAgent> neighbours) {
        Vektor force1 = new Vektor(0, 0, 0);
        Vektor force2 = new Vektor(0, 0, 0);
        ArrayList<Obstacle> visibleObstacles = new ArrayList<>();
        // Get visible obstacles
        for (Obstacle o : obstacles) {
            if (Vektor.subtract(agent.getPosition(), o.location).getLenght() - AGENT_RADIUS - o.radius < formulars.r_straight) {
                visibleObstacles.add(o);
            }
        }
        for (Obstacle o : visibleObstacles) {
            Vektor oPosition = formulars.sphericalObstaclePosition(o, agent);
            Vektor v = Vektor.subtract(agent.getPosition(), oPosition);
            Vektor f = formulars.n_ij(Vektor.subtract(agent.getPosition(), oPosition));
            f.multiply(formulars.Phi_beta(formulars.sigma_norm(v)));
            force1.add(f);
        }
        for (Obstacle o : visibleObstacles) {
            Vektor oPosition = formulars.sphericalObstaclePosition(o, agent);
            Vektor oVelocity = formulars.sphericalObstacleVelocity(o, agent);
            double d = formulars.b_ik(Vektor.subtract(agent.getPosition(), oPosition));
            Vektor v = Vektor.subtract(agent.getVelocity(), oVelocity);
            v.multiply(d);
            force2.add(v);
        }

        force1.multiply(formulars.c1_beta);
        force2.multiply(formulars.c2_beta);

        force1.add(force2);
        return force1;
    }

    /**
     * Gradient-based Term (desired distance)
     *
     * @param agent      current Agent
     * @param neighbours List with all Neighbours
     * @return Vektor force
     */
    private Vektor fig(IAgent agent, ArrayList<IAgent> neighbours) {
        Vektor force = new Vektor(0, 0, 0);
        for (IAgent neighbour : neighbours) {
            Vektor v = Vektor.subtract(agent.getPosition(), neighbour.getPosition());
            double value = formulars.sigma_norm(v);
            v = formulars.n_ij(v);
            v.multiply(formulars.Phi_alpha(value));
            force.add(v);
        }
        force.multiply(formulars.c1_alpha);
        return force;
    }

    /**
     * Consenus Term (velocity)
     *
     * @param agent      current Agent
     * @param neighbours List with all Neighbours
     * @return Vektor force
     */
    private Vektor fid(IAgent agent, ArrayList<IAgent> neighbours) {
        Vektor force = new Vektor(0, 0, 0);
        for (IAgent neighbour : neighbours) {
            Vektor v = Vektor.subtract(agent.getVelocity(), neighbour.getVelocity());
            v.multiply(formulars.a_ij(Vektor.subtract(agent.getPosition(), neighbour.getPosition())));
            force.add(v);
        }
        force.multiply(formulars.c2_alpha);
        return force;
    }


    private ArrayList<IAgent> getNeighboursFrom(IAgent agent) {
        ArrayList<IAgent> neighbourList = new ArrayList<>();
        for (IAgent other : agents) {
            double distance = Vektor.subtract(agent.getPosition(), other.getPosition()).getLenght();
            distance = distance - 2*AGENT_RADIUS;
            if (other != agent && distance <= formulars.r) {
                neighbourList.add(other);
            }
        }
        return neighbourList;
    }

    private void analyzeMinDistObst() {
        double min = HIGH_VALUE;
        for (Obstacle o : obstacles) {
            for (IAgent q : agents) {
                if (!q.isSelected()) {
                    double d = Vektor.subtract(o.location, q.getPosition()).getLenght();
                    d = d - o.radius - AGENT_RADIUS;
                    if (d < min) {
                        min = d;
                    }
                }
            }
        }
        if (min != HIGH_VALUE) {
            floatPublisherminDistObstacle.sendValue(min);
        }
    }

    private void analyzeMin() {
        double d = HIGH_VALUE;
        double h = HIGH_VALUE;
        for (IAgent q : agents) {
            for (IAgent other : agents) {
                if (q != other) {
                    Vektor v = Vektor.subtract(q.getPosition(), other.getPosition());
                    double dist = v.getLenght();
                    dist -= 2 * AGENT_RADIUS;         // beachte Quad Radius / Breite
                    if (dist < d) {
                        d = dist;
                    }
                }
            }
            if (q.getPosition().getZ() < h) h = q.getPosition().getZ();
        }

        if (d != HIGH_VALUE) {
            //System.out.println(d);
            floatPublisherminDist.sendValue(d);
        }
        if (h != HIGH_VALUE) {
            floatPublisherminZ.sendValue(h);
        }
    }

    private void analyzeAvgLeader() {
        if (leader != null) {
            int count = 0;
            double d = 0;
            for (IAgent q : agents) {
                if (q != leader) {
                    Vektor v = Vektor.subtract(q.getPosition(), leader.getPosition());
                    d += v.getLenght();
                    d -= 2 * AGENT_RADIUS;            // beachte Quad Radius / Breite
                    count++;
                }
            }

            if (count > 0) {
                d = d / count;
                floatPublisherAvgLeader.sendValue(d);
            }
        }

    }

    public void showParameterDialog() {
        new OlfatiSaberParameterDialogAlgorithm3(formulars);
    }


}
