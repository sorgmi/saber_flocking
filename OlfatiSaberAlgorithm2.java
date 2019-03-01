package de.dlr.swex.algorithm.saber;

import java.util.ArrayList;

import de.dlr.swex.Vektor;
import de.dlr.swex.algorithm.IAgent;
import de.dlr.swex.algorithm.parameter.OlfatiSaberParameterDialog;
import de.dlr.swex.ros.FloatPublisher;
import de.dlr.swex.ros.RosLoggerPublisher;
import rosgraph_msgs.Log;

/**
 * Olfati-Saber Flocking - Algorithm 2 without Obstacle Avoidance)
 *
 * @see <a href="http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.121.7027&rep=rep1&type=pdf">http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.121.7027&rep=rep1&type=pdf</a>
 */

public class OlfatiSaberAlgorithm2 {


    private ArrayList<IAgent> agents = new ArrayList<>();
    private IAgent leader;
    private OlfatiSaberFormulars formulars = new OlfatiSaberFormulars();

    private float HIGH_VALUE = 1000;
    private float AGENT_RADIUS = 0.25f;
    private FloatPublisher floatPublisherminDist = new FloatPublisher("swarmTool/debug/minDist");
    private FloatPublisher floatPublisherminZ = new FloatPublisher("swarmTool/debug/minZ");
    private FloatPublisher floatPublisherAvgLeader = new FloatPublisher("swarmTool/debug/avgDistLeader");
    private FloatPublisher floatPublisherDeviationEnergy = new FloatPublisher("swarmTool/debug/DeviationEnergy");



    public OlfatiSaberAlgorithm2(ArrayList<IAgent> list) {
        this.agents = list;
    }

    public void iterate() {
        //System.out.println(formulars.d_alpha);
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
            Vektor fiy = fiy(agent);

            //fig.multiply(1.75);

            acceleration.add(fig);
            acceleration.add(fid);
            acceleration.add(fiy);
            velocity.add(acceleration);
            velocity.limit(1.0f);
            Vektor pos = agent.getPosition();
            pos.add(velocity);
            pos.add(velocity);
            safetyBox(agent,pos);
            acceleration.multiply(0);
        }
        analyzeMin();
        analyzeAvgLeader();
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
     * Gradient-based Term (desired distance)
     * @param agent current Agent
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
        return force;
    }

    /**
     * Consenus Term (velocity)
     * @param agent current Agent
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
        return force;
    }

    /**
     * Group objective / Navigational feedback --> Follow Leader
     * @param agent current Agent
     * @return Vektor force
     */
    private Vektor fiy(IAgent agent){
       /* Vektor q = agent.getPosition();
        Vektor p = agent.getVelocity();
        q.multiply(-formulars.c1);
        p.multiply(formulars.c2);
        return Vektor.subtract(p,q); */

        if(leader != null){
            Vektor q = Vektor.subtract(leader.getPosition(),agent.getPosition());
            Vektor p = Vektor.subtract(leader.getVelocity(),agent.getVelocity());
            q.multiply(-formulars.c1);
            p.multiply(formulars.c2);
            q.sub(p);
            return q;
        }
        else {
            return new Vektor(0,0,0);
        }

    }

    private ArrayList<IAgent> getNeighboursFrom(IAgent agent) {
        ArrayList<IAgent> neighbourList = new ArrayList<>();
        for (IAgent other : agents) {
            double distance = Vektor.subtract(agent.getPosition(), other.getPosition()).getLenght(); // todo soll Agent-Radius beachtet werden?
            if (other != agent && distance <= formulars.r) { // todo use radius  r or d_alpha ?
                neighbourList.add(other);
            }
        }
        return neighbourList;
    }

    private void analyzeDeviationEnergy() {

        double energy = 0;
        ArrayList<Edge> edges = new ArrayList<>();

            for (IAgent i : agents) {
                ArrayList<IAgent> neighbours = new ArrayList<>();
                neighbours = getNeighboursFrom(i);
                for(IAgent j : neighbours){
                    if(j != i){
                        addEdgeIfNotExists(edges, i,j);
                        double l = Vektor.subtract(i.getPosition(),j.getPosition()).getLenght();
                        l -= formulars.d;
                        l = l*l;
                        energy += l;
                    }
                }
            }

        if(edges.size()>0){
            double f = edges.size()+1;
            energy = energy / f;
            floatPublisherDeviationEnergy.sendValue(energy);
            //System.out.println(edges.size());
        }


    }

    private void addEdgeIfNotExists(ArrayList<Edge> edges, IAgent a, IAgent b){
        boolean found = false;
        for(Edge e : edges){
            if( (e.a == a || e.b ==a) && (e.b == b || e.b ==a)){
                found = true;
            }

        }
        if(!found){
            Edge newEdge = new Edge();
            newEdge.a = a;
            newEdge.b = b;
            edges.add(newEdge);
        }
    }

    private class Edge{
        IAgent a,b;
    }

    private void analyzeMin() {
        double d = HIGH_VALUE;
        double h = HIGH_VALUE;
        for (IAgent q : agents) {
            for (IAgent other : agents) {
                if (q != other) {
                    Vektor v = Vektor.subtract(q.getPosition(), other.getPosition());
                    double dist = v.getLenght();
                    dist -= 2*AGENT_RADIUS;         // beachte Quad Radius / Breite
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
                    d -= 2*AGENT_RADIUS;            // beachte Quad Radius / Breite
                    count++;
                }
            }

            if (count > 0) {
                d = d / count;
                floatPublisherAvgLeader.sendValue(d);
            }
        }

    }

    public void showParameterDialog(){
        new OlfatiSaberParameterDialog(formulars);
    }


}
