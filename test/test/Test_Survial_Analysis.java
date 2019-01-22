/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.logging.Level;
import java.util.logging.Logger;
import person.AbstractIndividualInterface;
import population.MSMPopulation;
import population.person.MultiSiteMultiStrainPersonInterface;
import population.person.RelationshipPerson_MSM;
import relationship.RelationshipMap;
import relationship.SingleRelationship;
import sim.SinglePopRunnable;
import util.FileZipper;

/**
 *
 * @author Bhui
 */
public class Test_Survial_Analysis {

    public static void main(String[] arg) throws IOException, InterruptedException, ClassNotFoundException, FileNotFoundException, ExecutionException {
        Test_Survial_Analysis tsa = new Test_Survial_Analysis();
        tsa.runAnalysis();
    }

    public void runAnalysis() throws FileNotFoundException, IOException, InterruptedException, ExecutionException {

        String dirName = "C:\\Users\\Bhui\\OneDrive - UNSW\\MSM_MulitSite\\Test\\Import_1_No_Treatment_R";
        int NUM_THREAD = 1;

        int snapLimit = 4;

        File exinctAt = new File(dirName, "snapCount.obj_newStrainExinctAtSnapshot.csv");

        BufferedReader reader = new BufferedReader(new FileReader(exinctAt));
        String line;
        int[] extinctSnap = new int[1000];

        while ((line = reader.readLine()) != null) {
            String[] ent = line.split(",");
            int pt = Integer.parseInt(ent[0]);
            extinctSnap[pt] = Integer.parseInt(ent[1]);
        }

        reader.close();

        Future<Integer>[] numSecondary = new Future[extinctSnap.length];

        ExecutorService executor = null;
        int inExe = 0;

        GetNumberSecondary[] runnable = new GetNumberSecondary[extinctSnap.length];

        for (int i = 0; i < runnable.length; i++) {

            File basePop = new File(dirName, "export_21780");
            basePop = new File(basePop, "pop_" + i + ".zip");

            if (executor == null) {
                executor = Executors.newFixedThreadPool(NUM_THREAD);
            }

            runnable[i] = new GetNumberSecondary(i, basePop);

            if (extinctSnap[i] >= snapLimit) {

                numSecondary[i] = executor.submit(runnable[i]);
                inExe++;
                if (inExe == NUM_THREAD) {
                    executor.shutdown();
                    if (!executor.awaitTermination(2, TimeUnit.DAYS)) {
                        System.err.println("Thread time-out!");
                    }
                }
                executor = null;
                inExe = 0;
            }

        }

        if (executor != null) {
            executor.shutdown();
            if (!executor.awaitTermination(2, TimeUnit.DAYS)) {
                System.err.println("Thread time-out!");
            }
        }

        File numSecondaryFile = new File(dirName, "TSA_numSecondary.csv");
        PrintWriter wri = new PrintWriter(numSecondaryFile);

        for (int i = 0; i < numSecondary.length; i++) {
            wri.print(i);
            wri.print(',');
            if(numSecondary[i] != null){
                wri.print(numSecondary[i].get());
            } else{
                wri.print(0);
            }
                   
            wri.println();
        }

        wri.close();

    }

    private class GetNumberSecondary implements Callable<Integer> {

        File popFile;
        int id;

        public GetNumberSecondary(int threadId, File popFile) {
            this.popFile = popFile;
            this.id = threadId;
        }

        @Override
        public Integer call() throws Exception {
            System.out.println("#" + id + ": Calculating # secondary case for " + popFile.getAbsolutePath());

            File temp = FileZipper.unzipFile(popFile, popFile.getParentFile());

            MSMPopulation pop;
            try (ObjectInputStream inStream = new ObjectInputStream(new FileInputStream(temp))) {
                pop = MSMPopulation.importMSMPopulation(inStream);
            }
            temp.delete();

            RelationshipPerson_MSM patient_zero = null;
            int secondary_case = 0;

            HashMap<Integer, RelationshipPerson_MSM> mapping = new HashMap<>();

            // Detecting patient zero
            if (pop != null) {

                pop.initialiseInfection(0); // 0 = using orginal RNG

                for (AbstractIndividualInterface person : pop.getPop()) {

                    RelationshipPerson_MSM p = (RelationshipPerson_MSM) person;

                    mapping.put(p.getId(), p);

                    boolean hasNewStrain = hasNewStrain(p);
                    if (hasNewStrain) {
                        patient_zero = p;
                        System.out.println("#" + id + ": Patient zero id = " + p.getId());
                    }
                }

                if (patient_zero != null) {

                    while (hasNewStrain(patient_zero)) // One sim run                
                    {

                        HashMap<Integer, int[]> partnersCollection = new HashMap<>();

                        for (RelationshipMap relMap : pop.getRelMap()) {
                            if (relMap.containsVertex(patient_zero.getId())) {
                                SingleRelationship[] relArr = relMap.edgesOf(patient_zero.getId()).toArray(new SingleRelationship[0]);
                                fillPartnerCollection(relArr, patient_zero, partnersCollection, mapping);
                            }
                        }

                        pop.advanceTimeStep(1);

                        for (RelationshipMap relMap : pop.getRelMap()) {
                            if (relMap.containsVertex(patient_zero.getId())) {
                                SingleRelationship[] relArr = relMap.edgesOf(patient_zero.getId()).toArray(new SingleRelationship[0]);
                                secondary_case += getNumSecondary(relArr, patient_zero, mapping, partnersCollection);
                            }
                        }

                    }

                }
            }

            return secondary_case;

        }

        protected int getNumSecondary(SingleRelationship[] relArr,
                RelationshipPerson_MSM patient_zero,
                HashMap<Integer, RelationshipPerson_MSM> mapping,
                HashMap<Integer, int[]> partnersCollection) {

            int newSecondary = 0;
            for (SingleRelationship rel1 : relArr) {
                int[] partners = rel1.getLinksValues();

                RelationshipPerson_MSM partner;

                if (patient_zero.getId() == partners[0]) {
                    partner = mapping.get(partners[1]);
                } else {
                    partner = mapping.get(partners[0]);
                }

                if (hasNewStrain(partner)) {
                    if (partnersCollection.containsKey(partner.getId())) {
                        boolean newlyGained = false;
                        int[] strainStat = partner.getCurrentStrainsAtSite();
                        int[] pastStat = partnersCollection.get(partner.getId());

                        for (int s = 0; s < pastStat.length; s++) {
                            newlyGained |= (strainStat[s] & 0b10) > 0 && (pastStat[s] & 0b10) == 0;
                        }

                        if (newlyGained) {
                            newSecondary++;
                        }
                    } else {
                        newSecondary++;
                    }

                }

            }

            return newSecondary;
        }

        protected void fillPartnerCollection(
                SingleRelationship[] relArr,
                RelationshipPerson_MSM patient_zero,
                HashMap<Integer, int[]> partnersCollection,
                HashMap<Integer, RelationshipPerson_MSM> mapping) {

            for (SingleRelationship rel1 : relArr) {
                int[] partners = rel1.getLinksValues();

                RelationshipPerson_MSM partner;

                if (patient_zero.getId() == partners[0]) {
                    partner = mapping.get(partners[1]);
                } else {
                    partner = mapping.get(partners[0]);
                }
                partnersCollection.put(partner.getId(), Arrays.copyOf(partner.getCurrentStrainsAtSite(),
                        partner.getCurrentStrainsAtSite().length));
            }
        }

        protected boolean hasNewStrain(RelationshipPerson_MSM p) {
            int[] strainStat = p.getCurrentStrainsAtSite();
            boolean hasNewStrain = false;
            for (int i = 0; i < strainStat.length; i++) {
                hasNewStrain |= (strainStat[i] & 0b10) > 0;
            }
            return hasNewStrain;
        }

    }
}
