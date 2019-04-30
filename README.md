### GALATI Mathias 
### Etudiant à l'Unisversité de Strasbourg en Bioligie Structurale intégrative et Bio-Informatique
### Stage M2S3 2019 : CIRAD Montpellier
***************
### Préambule
Ce dépot GitHub permet de regrouper l'ensemble des démarches et protocoles adoptées danns le cadre d'un sujet de stage en métabarcoding et d'un autre sujet en métagénomique.
***************
### Sommaire :
1. Descritption rapide de l'ensemble des répertoires du GitHub ainsi que l'architecture adoptée
2. Guide utilisateur des scripts correspondants aux workflows Métabarcoding
3. Guide utilisateur des scripts correspondants aux workflows Métagénomique
4. Ensemble des ressources utilisées
<!-- -->
Un fichier de biliographie sera édité et mis à dispostion en fin de stage
**************
### 1. Architucture des répertoires :
- Biblio : Fichier qui contiendra la bilbiographie de l'ensemble du stage
- metabarcoding/           
    * Répertoire contenant les différents workflows de scripts permettant l'analyse des données 
- metabarcoding/16S/    
    * Répertoire contentenant les différents scripts pour données taxonomique de type 16S :  
        * 0 = Sous-échantillonnage permettant de dévelloper les scripts avec des temps de calculs courts
        * 1 = Analyse fastQC / multiQC
        * 2 = Suppression des primers et adaptateurs des séquences brut de séquencage avec Trim Galore et Cutadapt
        * 3.1 = Analyse avec l'approche dada2
        * 3.2 = Analyse avec l'approche deblur
        * 3.3 = Analyse avec l'approche vsearch
        
- metabarcoding/ITS/    
    * Répertoire contentenant les différents scripts pour des données taxonomique de type ITS :  
        * 0 = Sous-échantillonnage permettant de dévelloper les scripts avec des temps de calculs courts
        * 1 = Analyse fastQC / multiQC
        * 2 = Suppression des primers et adaptateurs des séquences brut de séquencage avec Trim Galore et Cutadapt
        * 3.1 = Analyse avec l'approche dada2
        * 3.2 = Analyse avec l'approche deblur
        * 3.3 = Analyse avec l'approche vsearch  
- metagenomique/    
    * Répertoire contentenant :  
        * Un fichier de texte de notes sur le début d'élaboration d'un workflow avec Snakemake sous conda 

- tuto
    * Répertoire contentenant :  
        * Un fichier de texte de notes sur le début d'élaboration d'un workflow avec Snakemake sous conda 
*************** 
### 2. Métabarcoding    

***************
### 3. Métagénomique    
***************
### 4. Ressources    
[WORK IN PROGRESS]
