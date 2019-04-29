package org.noise_planet.nm_run_multiples_instances

/** author : Aumond Pierre

 Ce code permet de faire une unique de simulation

 Le résultat se trouve dans un dossier zippé
 Il s'agit d'une version Postgre/Postgis (lecture/ecriture plus lent qu'avec h2gis mais jointures beaucoup plus rapides)

 **/

// Importation des librairies nécessaire au code


/** author : Aumond Pierre

 Ce code permet de faire une unique de simulation

 Le résultat se trouve dans un dossier zippé
 Il s'agit d'une version Postgre/Postgis (lecture/ecriture plus lent qu'avec h2gis mais jointures beaucoup plus rapides)

 **/

// Importation des librairies nécessaire au code
import org.noise_planet.noisemodelling.emission.EvaluateRoadSourceCnossos
import org.noise_planet.noisemodelling.emission.RSParametersCnossos
import org.noise_planet.noisemodelling.propagation.ComputeRays

class OnlyEmission {
    static void main(String[] args) {
        OnlyEmission oneRunEmission = new OnlyEmission()
        oneRunEmission.run()
    }

    void run() {
/**
 ///////////////////////////////////////////
 // Paramètres d'entrée et initialisations //
 ///////////////////////////////////////////

 */
        String workspace_output = "D:\\aumond\\Documents\\Recherche_hors_projet\\2019_04_PasseportRecherche\\Traitement"

        String input_file = "D:\\aumond\\Documents\\Recherche_hors_projet\\2019_04_PasseportRecherche\\Donnees\\Eleves6.csv"

        /**
         ///////////////////////////////////////////
         // FIN Paramètres d'entrée et initialisations //
         ///////////////////////////////////////////

         */


        try {
            PrintWriter pw = new PrintWriter(new File("D:\\sortieModelFlow2.csv"))
            System.out.println("Init")

            /*VL;ML;PL;Vit;Vitest;NbVoies;Dist;Haute*/

            int nvar = 8

            ArrayList<Double> VL = new ArrayList<Double>()
            ArrayList<Double> ML = new ArrayList<Double>()
            ArrayList<Double> PL = new ArrayList<Double>()
            ArrayList<Double> Vit = new ArrayList<Double>()
            ArrayList<Double> Vitest = new ArrayList<Double>()
            ArrayList<Double> NbVoies = new ArrayList<Double>()
            ArrayList<Double> Dist = new ArrayList<Double>()
            ArrayList<Double> Haute = new ArrayList<Double>()
            ArrayList<Double> Temps = new ArrayList<Double>()
            ArrayList<Double> Vlh = new ArrayList<Double>()
            ArrayList<Double> MlH = new ArrayList<Double>()
            ArrayList<Double> Plh = new ArrayList<Double>()
            ArrayList<Double> PK = new ArrayList<Double>()
            ArrayList<Double> Result = new ArrayList<Double>()
            def list = [63, 125, 250, 500, 1000, 2000, 4000, 8000]


            new File(input_file).splitEachLine(";") { fields ->
                def res = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                def kk = 0
                VL.add(fields[0].toFloat())
                ML.add(fields[1].toFloat())
                PL.add(fields[2].toFloat())
                Vit.add(fields[3].toFloat())
                Vitest.add(fields[4].toFloat())
                NbVoies.add(fields[5].toFloat())
                Dist.add(fields[6].toFloat())
                Haute.add(fields[7].toFloat())
                Temps.add(fields[8].toFloat())
                Vlh.add(fields[9].toFloat())
                MlH.add(fields[10].toFloat())
                Plh.add(fields[11].toFloat())
                PK.add(fields[12].toFloat())

                for (f in list) {
                    //RSParametersCnossos srcParameters = new RSParametersCnossos(fields[4].toFloat(), fields[4].toFloat(), fields[4].toFloat(), 0, 0,
                    //        fields[0].toFloat(), fields[1].toFloat(), fields[2].toFloat(), 0, 0, f, 20, "NL01", 0, 0, 250, 1)
                    //RSParametersCnossos srcParameters = new RSParametersCnossos(fields[4].toFloat(), fields[4].toFloat(), fields[4].toFloat(), 0, 0,
                    //        3600, 0, 0, 0, 0, f, 20, "NL01", 0, 0, 250, 1)
                    RSParametersCnossos srcParameters = new RSParametersCnossos(fields[4].toFloat(), fields[4].toFloat(), fields[4].toFloat(), 0, 0,
                            fields[9].toFloat(), fields[10].toFloat(), fields[11].toFloat(), 0, 0, f, 20, "NL01", 0, 0, 250, 1)


                    srcParameters.setSlopePercentage(RSParametersCnossos.computeSlope(0, 0, 10))
                    res[kk] = EvaluateRoadSourceCnossos.evaluate(srcParameters)
                    res[kk] = 10 * Math.log10(Math.pow(10, EvaluateRoadSourceCnossos.evaluate(srcParameters) / 10) + Math.pow(10, res[kk] / 10))
                    kk++

                }

                res[8] = ComputeRays.wToDba(ComputeRays.dbaToW(res[0] - 26.2) + ComputeRays.dbaToW(res[1] - 16.1) + ComputeRays.dbaToW(res[2] - 8.6) + ComputeRays.dbaToW(res[3] - 3.2) + ComputeRays.dbaToW(res[4]) + ComputeRays.dbaToW(res[5] + 1.2) + ComputeRays.dbaToW(res[6] + 1) + ComputeRays.dbaToW(res[7] - 1.1))


                Result.add(res)
                pw.write(fields[12].toFloat() + ";" + res[8].toString() + "\n")

            }

            pw.close()


        } finally {

        }
    }
}