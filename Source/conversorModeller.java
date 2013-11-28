/*
 * Script para converter sequencias para um formato no qual o modeller consegue
 * utilizar para modelar proteínas.
 */
package tcc;

import java.io.File;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * @author Goron
 */
public class conversorModeller {

    /**
     * Salva o conteúdo de uma variável em um arquivo
     * @param file
     * @param content
     * @param add se true adicionar no final do arquivo
     * @throws IOException 
     */
    public static void save(String file, String content, boolean add) throws IOException {

        // Abre um arquivo para escrita verificando se será para sobrescrever ou não
        FileWriter fw = new FileWriter(file, add);

        // Escreve no arquivo
        fw.write(content);

        //Fecha o arquivo
        fw.close();
    }
    // Variável global com informações para serem postas no arquivo
    static String[][] splited2;
    static List<String> listSeq = new ArrayList<String>();

    /**
     * Faz todo o trabalho sobre o arquivo .pdb chamando as outras funções para
     * serem executadas, será a função chamada pelo botão quando houver
     * interface gráfica
     */
    public static void workerPdb(String Path) {

        // String para armazenar o texto que retorna da função load
        String fileText;

        // Listas de string para armazenar o conteudo dos arquivos apos o split
        String[] splited;

        //Variável para armazenar a quantidade de cadeias
        int ChainCount = 0;

        try {



            // Carregar o PDB
            fileText = load(Path);

            // Faz o split no PDB reutilizando a variável splited           
            splited = fileText.split("ATOM      1");
            splited[1] = "ATOM      1" + splited[1];
            splited = splited[1].split("TER");
            String[] auxiliar = {""};
            ChainCount = splited.length - 2; // -2 pois descarta a posição 0, que é vazia e a última que não possui informações úteis
            splited2 = new String[splited.length][3];

            for (int i = 0; i < ChainCount; i++) {
                auxiliar = splited[i].split("ATOM");
                splited2[i][0] = auxiliar[1].split("\\s+")[4];
                splited2[i][1] = auxiliar[1].split("\\s+")[5];
                splited2[i][2] = auxiliar[auxiliar.length - 2].split("\\s+")[5];
                System.out.println(splited2[i][0] + " " + splited2[i][1] + " " + splited2[i][2] + " ");
            }

        } catch (FileNotFoundException ex) {

            // Exceção de arquivo não existente
            Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, ex);

        } catch (IOException ex) {

            // Exceção de entrada e saída
            Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, ex);

        }
    }

    /**
     * Faz todo o trabalho sobre o arquivo .aa chamando as outras funções para
     * serem executadas, será a função chamada pelo botão quando houver
     * interface gráfica
     */
    public static List<String> workerAa(String path, int times) {

        // Tokenizers para fazer split
        StringTokenizer text = null;
        StringTokenizer splitedText = null;

        // String para armazenar o texto que retorna da função load
        String fileText = "";

        // Listas de string para armazenar o conteudo dos arquivos apos o split
        // e para escrever a lista de retorno que vai ser mostrada no jList
        List<String> list = new ArrayList<String>();
        List<String> seqNames = new ArrayList<String>();
        List<String> auxList = new ArrayList<String>();

        try {

            // Limpa as listas
            list.clear();
            seqNames.clear();
            
            // Carregar a sequencia .aa
            fileText = load(path);

            // Faz o split no texto .aa
            text = Splitter(fileText, ">");

            // Armazena as partes da sequencia em uma lista
            while (text.hasMoreTokens()) {
                list.add(text.nextToken());
            }

            // Armazena a lista em outra lista apos fazer split novamente
            for (int i = 0; i < list.size(); i++) {
                splitedText = Splitter(list.get(i), " . ");
                while (splitedText.hasMoreTokens()) {
                    listSeq.add(splitedText.nextToken());
                }
            }

            for (int i = 0; i < listSeq.size(); i = i + 2) {
                auxList.add(listSeq.get(i).toString());
            }

            for(int i=0;i<auxList.size()/times;i++){
                seqNames.add(auxList.get(i));
            }
            
            return seqNames;
        } catch (FileNotFoundException ex) {

            // Exceção de arquivo não existente
            Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, ex);

        } catch (IOException ex) {

            // Exceção de entrada e saída
            Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, ex);

        }
        return null;
    }

    public static void Start(String mySequence, boolean selected, int value) throws IOException {

        // Preciso contar decentemente os residuos da minha proteína, e adicionar suporte pra caso eu tenha
        // juntado todas as sequencias (aquele erro comum que temos de ter que juntar todas sequencias do PDB)
        // as sequencias juntas vão dentro do "if(selected)" e a contagem de resíduos vai em residues.length()
        
        String output = "";
        String residues;
        int residueCount = 0;
        int notvalue = 0;
        
        if(value==0)
        {
            notvalue = 1;
        }
        
        residues = listSeq.get(1).replaceAll("-", "");
        residues = residues.replaceAll("\\r?\\n","");

        if (selected) {
            
            int tamanho;
            tamanho = 0;
            
            for(int i = 0; i<splited2.length-2; i++)
            {
                tamanho += Integer.parseInt(splited2[i][2]) - Integer.parseInt(splited2[i][1]);
            }
            
           tamanho += Integer.parseInt(splited2[0][1]);
           
           if(value==1)
           {
            output = ">P1;" + listSeq.get(notvalue*2)
                    + "\r\nstructureX:" + listSeq.get(notvalue*2)
                    + ": " + splited2[0][1] + "     :" + splited2[0][0] + ":" + tamanho + "  :" + splited2[0][0] + "::: :\r\n" + listSeq.get(1)
                    + "*\r\n\r\n>P1;" + listSeq.get(value*2)
                    + "\r\nsequence:" + listSeq.get(value*2)
                    + ":1    :A:" + residues.length() + "   :A::: :\r\n" + listSeq.get(3) + "*";
           } else {
               
                    output = ">P1;" + listSeq.get(notvalue*2)
                    + "\r\nstructureX:" + listSeq.get(notvalue*2)
                    + ": " + splited2[0][1] + "     :" + splited2[0][0] + ":" + tamanho + "  :" + splited2[0][0] + "::: :\r\n" + listSeq.get(3)
                    + "*\r\n\r\n>P1;" + listSeq.get(value*2)
                    + "\r\nsequence:" + listSeq.get(value*2)
                    + ":1    :A:" + residues.length() + "   :A::: :\r\n" + listSeq.get(1) + "*";
              
           }
           
        } else {
            
            if(value == 1){
            output = ">P1;" + listSeq.get(notvalue*2)
                    + "\r\nstructureX:" + listSeq.get(notvalue*2)
                    + ": " + splited2[0][1] + "     :" + splited2[0][0] + ":" + splited2[0][2] + "  :" + splited2[0][0] + "::: :\r\n" + listSeq.get(1)
                    + "*\r\n\r\n>P1;" + listSeq.get(value*2)
                    + "\r\nsequence:" + listSeq.get(value*2)
                    + ":1    :A:" + residues.length() + "   :A::: :\r\n" + listSeq.get(3) + "*";
            } else {
                
                output = ">P1;" + listSeq.get(notvalue*2)
                    + "\r\nstructureX:" + listSeq.get(notvalue*2)
                    + ": " + splited2[0][1] + "     :" + splited2[0][0] + ":" + splited2[0][2] + "  :" + splited2[0][0] + "::: :\r\n" + listSeq.get(3)
                    + "*\r\n\r\n>P1;" + listSeq.get(value*2)
                    + "\r\nsequence:" + listSeq.get(value*2)
                    + ":1    :A:" + residues.length() + "   :A::: :\r\n" + listSeq.get(1) + "*";
            }
        }
        // Print de teste
        System.out.println(output);

        // Salva o arquivo de texto
        save("align.ali", output, false);
    }
    
     public static void StartModel(String mySequence, boolean selected, int value, int templates) throws IOException {

               
        String output = "";
        
        int notvalue = 0;
        
        if(value==0)
        {
            notvalue = 1;
        }
      

           output = "# Homology modeling with multiple templates\n" +
                     "from modeller import *              # Load standard Modeller classes\n" +
                     "from modeller.automodel import *    # Load the automodel class\n" +

                     "log.verbose()    # request verbose output" +
                     "env = environ()  # create a new MODELLER environment to build this model in\n" +

                     "# directories for input atom files\n" +
                     "env.io.atom_files_directory = ['.', '../atom_files']\n" +

                     "a = automodel(env,\n" +
                     "              alnfile  = 'align-multiple.ali', # alignment\n filename\n" +
                     "              knowns   = ('"+ listSeq.get(notvalue*2).substring(0, 4) +"'),     # codes of\n the templates\n" +
                     "              sequence = '" + listSeq.get(value*2) + "')               # code of the\n target\n" +
                     "a.starting_model= 1                 # index of the first model\n" +
                     "a.ending_model  = " + templates + "                 # index of the last model\n" +
                                    "# (determines how many models to calculate)\n" +
                     "a.make()                            # do the actual homology modeling\n";
           

            // Print de teste
        System.out.println(output);

        // Salva o arquivo de texto
        save("model.py", output, false);
    }

    /**
     * Carrega o conteúdo de um arquivo em uma stirng, se o arquivo não existir,
     * retornará null
     * @param file
     * @return butOut
     * @throws FileNotFoundException
     * @throws IOException 
     */
    public static String load(String file) throws FileNotFoundException, IOException {

        // Abre um arquivo
        File opened = new File(file);

        // Checa se o arquivo existe, se não existir retorna nulo
        if (!opened.exists()) {
            return null;
        }

        // Cria um buffer de leitura
        BufferedReader br = new BufferedReader(new FileReader(file));

        // Cria uma string builder para concatenar as linhas
        StringBuilder bufOut = new StringBuilder();

        // Cria uma string para ler a linha
        String line;

        // Concatena as strings no string builder
        while ((line = br.readLine()) != null) {
            bufOut.append(line).append("\n");
        }

        // Fecha o arquivo
        br.close();

        // Retorna o string builderm o transformando em uma string
        return bufOut.toString();
    }

    /**
     * Faz um split usando o tokenizer
     * @param text
     * @param splitter
     * @return tokens
     */
    public static StringTokenizer Splitter(String text, String splitter)
            throws FileNotFoundException, IOException {
        StringTokenizer tokens = new StringTokenizer(text, splitter);
        return tokens;
    }
}
