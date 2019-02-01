/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

/**
 *
 * @author bbcbh
 */
public class Test_SubAnalysis_Abstract {

    public static void main(String[] arg) throws IOException {
        File baseDir = new File("C:\\Users\\bbcbh\\Downloads\\B2");

        File newStrainExtinctCSV = new File(baseDir, "snapCount.obj_newStrainExinctAtSnapshot.csv");
        File popSnapDir = new File(baseDir, "export_7920");
        File newStrainDir = new File(baseDir, "newStrainSpread");

        int[] PERSIST_RANGE = new int[]{2, 4};
        int[] GROUPING_CASUSAL_PARNTERS = new int[]{0, 4, 9, 19, 39};

        int[][] matrix_persist = readCSVFile(newStrainExtinctCSV, 0);
        int numSimTotal = popSnapDir.listFiles(new FileFilter() {
            @Override
            public boolean accept(File pathname) {
                return pathname.getName().endsWith(".csv");
            }
        }).length;

        numSimTotal = Math.min(numSimTotal, matrix_persist.length - 1);

        int[][][] index_stat = new int[PERSIST_RANGE.length][GROUPING_CASUSAL_PARNTERS.length + 1][2]; // yes, total
        int[][][] index_partner_stat = new int[PERSIST_RANGE.length][GROUPING_CASUSAL_PARNTERS.length + 1][2]; // yes, total

        for (int sN = 0; sN < numSimTotal; sN++) {
            File snapCSV = new File(popSnapDir, "IndivdualSnap_" + sN + ".csv");

            System.out.println("Analysing " + snapCSV.getAbsolutePath());
            int[][] matrix_popSnap = readCSVFile(snapCSV, 1);

            int[] persistUpTo = matrix_persist[sN];

            if (sN != persistUpTo[0]) {
                System.err.println("Row mismatch at Row #" + sN + ". Attemping full search");

                for (int[] p : matrix_persist) {
                    if (p[0] == sN) {
                        persistUpTo = p;
                        break;
                    }
                }
            }

            int[] patient_zero = null;
            HashMap<Integer, int[]> regPartnerRec = new HashMap();

            File strainSpread = new File(newStrainDir, "newStrainSpread_" + sN + "_7919.csv");

            if (strainSpread.exists()) {
                BufferedReader r = new BufferedReader(new FileReader(strainSpread));
                String line;

                boolean partnerEntNext = false;
                while ((line = r.readLine()) != null) {
                    if (partnerEntNext) {
                        do {

                            String[] ent = line.split(",");
                            if (ent[3].equals("0")) {
                                regPartnerRec.put(Integer.parseInt(ent[0]), new int[0]);
                            }
                        } while ((line = r.readLine()) != null && !line.startsWith("Time"));

                    }
                    partnerEntNext = line != null && line.startsWith("PATIENT_ZERO_PARNTER_ID");

                }
            }

            for (int[] popEnt : matrix_popSnap) {
                if (popEnt[11] == 2) {
                    patient_zero = Arrays.copyOf(popEnt, popEnt.length);
                }

                if (regPartnerRec.containsKey(popEnt[0])) {
                    regPartnerRec.put(popEnt[0], popEnt);
                }

            }

            if (patient_zero != null) {

                int cPt;

                cPt = Arrays.binarySearch(GROUPING_CASUSAL_PARNTERS, patient_zero[5]);

                if (cPt < 0) {
                    cPt = -(cPt + 1);
                }

                for (int rangePt = 0; rangePt < PERSIST_RANGE.length; rangePt++) {
                    index_stat[rangePt][cPt][1]++;
                    if (persistUpTo[1] > PERSIST_RANGE[rangePt]) {
                        index_stat[rangePt][cPt][0]++;
                    }
                }

                for (int[] partner_info : regPartnerRec.values()) {
                    if (partner_info.length != 0) {
                        cPt = Arrays.binarySearch(GROUPING_CASUSAL_PARNTERS, partner_info[5]);

                        if (cPt < 0) {
                            cPt = -(cPt + 1);
                        }

                        for (int rangePt = 0; rangePt < PERSIST_RANGE.length; rangePt++) {
                            index_partner_stat[rangePt][cPt][1]++;
                            if (persistUpTo[1] > PERSIST_RANGE[rangePt]) {
                                index_partner_stat[rangePt][cPt][0]++;
                            }
                        }
                    }

                }

            }
        }

        for (int rangePt = 0; rangePt < PERSIST_RANGE.length; rangePt++) {

            int[][] statEnt = index_stat[rangePt];
            System.out.println("Persist > " + (PERSIST_RANGE[rangePt] / 4.0f));

            for (int cPt = 0; cPt < statEnt.length; cPt++) {
                if (cPt < GROUPING_CASUSAL_PARNTERS.length) {
                    System.out.println("<=" + (GROUPING_CASUSAL_PARNTERS[cPt]));
                } else {
                    System.out.println(">=" + (GROUPING_CASUSAL_PARNTERS[GROUPING_CASUSAL_PARNTERS.length - 1] + 1));
                }

                float roundPercentage = Math.round((1000f * statEnt[cPt][0]) / statEnt[cPt][1]) / 10f;

                System.out.println("     " + statEnt[cPt][0] + " / " + statEnt[cPt][1]
                        + " (" + roundPercentage + "%)");
            }

            System.out.println("Reg stat");
            System.out.println("Persist > " + (PERSIST_RANGE[rangePt] / 4.0f));

            statEnt = index_partner_stat[rangePt];
            for (int cPt = 0; cPt < statEnt.length; cPt++) {
                if (cPt < GROUPING_CASUSAL_PARNTERS.length) {
                    System.out.println("<=" + (GROUPING_CASUSAL_PARNTERS[cPt]));
                } else {
                    System.out.println(">=" + (GROUPING_CASUSAL_PARNTERS[GROUPING_CASUSAL_PARNTERS.length - 1] + 1));
                }

                float roundPercentage = Math.round((1000f * statEnt[cPt][0]) / statEnt[cPt][1]) / 10f;

                System.out.println("     " + statEnt[cPt][0] + " / " + statEnt[cPt][1]
                        + " (" + roundPercentage + "%)");
            }

        }

    }

    private static int[][] readCSVFile(File csv, int skip) throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader(csv));
        String line;

        ArrayList<int[]> rows = new ArrayList();
        // Skip
        for (int i = 0; i < skip; i++) {
            line = reader.readLine();
        }

        // First line;
        while ((line = reader.readLine()) != null) {
            String[] ent = line.split(",");
            int[] rEntry = new int[ent.length];
            for (int i = 0; i < rEntry.length; i++) {
                rEntry[i] = Integer.parseInt(ent[i]);
            }
            rows.add(rEntry);
        }

        int[][] res = new int[rows.size()][];
        int r = 0;

        for (int[] rEntry : rows) {
            res[r] = rEntry;
            r++;
        }
        return res;

    }

}
