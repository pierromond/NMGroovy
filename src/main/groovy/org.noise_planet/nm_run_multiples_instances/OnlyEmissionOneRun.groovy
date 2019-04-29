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
import org.noise_planet.noisemodelling.emission.EvaluateRoadSourceDynamic
import org.noise_planet.noisemodelling.emission.RSParametersCnossos
import org.noise_planet.noisemodelling.emission.RSParametersDynamic
import org.noise_planet.noisemodelling.propagation.ComputeRays

class OnlyEmissionOneRun {
    static void main(String[] args) {
        OnlyEmissionOneRun oneRunEmission = new OnlyEmissionOneRun()
        oneRunEmission.run()
    }

    void run() {
/**
 ///////////////////////////////////////////
 // Paramètres d'entrée et initialisations //
 ///////////////////////////////////////////

 */
        /**
         ///////////////////////////////////////////
         // FIN Paramètres d'entrée et initialisations //
         ///////////////////////////////////////////

         */


        try {
            PrintWriter pw = new PrintWriter(new File("D:\\sortieModelFlow2.csv"))
            System.out.println("Init")

            /*VL;ML;PL;Vit;Vitest;NbVoies;Dist;Haute*/

            int kk = 0
            ArrayList<Double> res = new ArrayList<>()
            ArrayList<Double> res2 = new ArrayList<>()
            def list = [63, 125, 250, 500, 1000, 2000, 4000, 8000]

            for (f in list) {
                RSParametersDynamic srcParameters = new RSParametersDynamic(60, 0, 0, 0, f, 20, 1, false, 2000, 0,0, 1)
                res[kk] = EvaluateRoadSourceDynamic.evaluate(srcParameters)

                RSParametersCnossos srcParameters2 = new RSParametersCnossos(60, 0, 0, 0, 0,
                       3600, 0, 0, 0, 0, f, 20, "NL01", 0, 0, 250, 1)

                res2[kk] = 10 * Math.log10(Math.pow(10, EvaluateRoadSourceCnossos.evaluate(srcParameters2) / 10) + Math.pow(10, res[kk] / 10))
                kk++

                }

            res[8] = ComputeRays.wToDba(ComputeRays.dbaToW(res[0] ) + ComputeRays.dbaToW(res[1] ) + ComputeRays.dbaToW(res[2] ) + ComputeRays.dbaToW(res[3] ) + ComputeRays.dbaToW(res[4]) + ComputeRays.dbaToW(res[5] ) + ComputeRays.dbaToW(res[6] ) + ComputeRays.dbaToW(res[7] ))
            res2[8] = ComputeRays.wToDba(ComputeRays.dbaToW(res[0] ) + ComputeRays.dbaToW(res[1] ) + ComputeRays.dbaToW(res[2] ) + ComputeRays.dbaToW(res[3] ) + ComputeRays.dbaToW(res[4]) + ComputeRays.dbaToW(res[5] ) + ComputeRays.dbaToW(res[6] ) + ComputeRays.dbaToW(res[7] ))

            res[9] = ComputeRays.wToDba(ComputeRays.dbaToW(res[0] - 26.2) + ComputeRays.dbaToW(res[1] - 16.1) + ComputeRays.dbaToW(res[2] - 8.6) + ComputeRays.dbaToW(res[3] - 3.2) + ComputeRays.dbaToW(res[4]) + ComputeRays.dbaToW(res[5] + 1.2) + ComputeRays.dbaToW(res[6] + 1) + ComputeRays.dbaToW(res[7] - 1.1))
            res2[9] = ComputeRays.wToDba(ComputeRays.dbaToW(res[0] - 26.2) + ComputeRays.dbaToW(res[1] - 16.1) + ComputeRays.dbaToW(res[2] - 8.6) + ComputeRays.dbaToW(res[3] - 3.2) + ComputeRays.dbaToW(res[4]) + ComputeRays.dbaToW(res[5] + 1.2) + ComputeRays.dbaToW(res[6] + 1) + ComputeRays.dbaToW(res[7] - 1.1))


            pw.write(res[8].toString() + "\n")
            pw.write(res2[8].toString() + "\n")
             pw.close()


        } finally {

        }
    }
}