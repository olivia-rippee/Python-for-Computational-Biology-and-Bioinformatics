import os
os.chdir("C:/Users/olivi/OneDrive/Python/Bioinformatics IV - Molecular Evolution/Data")


# -----------------------------------------------
# Expected Number of Peptide Matches
# -----------------------------------------------

def ExpectedPeptideMatches(probability, proteome_length):
    '''Calculate the expected number of peptide matches in a decoy proteome.
    E = N * p
    
    Input: Probability of a peptide matching (0 < p < 1) and length of the decoy proteome.
    Output: Expected number of peptide matches.
    '''
    expected_matches = probability * proteome_length
    return expected_matches


# Say that the probability of a collection of peptides Dictionary is equal to 
# 0.00012. What is the expected approximate number of peptide matches from 
# Dictionary that we expect to find in a decoy proteome of length 500,000?
# ------------------------------------------------------------------------
probability, proteome_length = 0.00012, 500000
print(ExpectedPeptideMatches(probability, proteome_length))  # Output: 60



# -----------------------------------------------
# Top Score Peptide from Proteome
# -----------------------------------------------

def PeptideIdentification(Spectrum, Proteome, MassTable):
    '''Find substring of Proteome with max score against SpectrumVector.'''
    spectrum = list(map(int, Spectrum.strip().split()))
    max_score = float('-inf')
    best_peptide = ""
    
    for start in range(len(Proteome)):
        mass_sum = 0
        for end in range(start, len(Proteome)):
            mass_sum += MassTable[Proteome[end]]
            if mass_sum > len(spectrum):
                break
            peptide = Proteome[start:end+1]
            s = Score(peptide, spectrum)
            if s > max_score:
                max_score = s
                best_peptide = peptide
    return best_peptide

def PeptideScoreVector(Peptide, MassTable):
    '''Convert peptide to its prefix-mass vector representation.'''
    prefix_mass = 0
    vector = []
    for aa in Peptide:
        prefix_mass += MassTable[aa]
        vector.append(prefix_mass)
    return vector

def Score(Peptide, Spectrum):
    '''Score a peptide against a spectral vector.'''
    vec = PeptideScoreVector(Peptide, MassTable)
    score = 0
    for mass in vec:
        if mass <= len(Spectrum):
            score += Spectrum[mass-1]  # 0-indexed
    return score



# Example 1
# -----------
Spectrum = "0 0 0 4 -2 -3 -1 -7 6 5 3 2 1 9 3 -8 0 3 1 2 1 8"
Proteome = "XZZXZXXXZXZZXZXXZ"
MassTable = {"X":4, "Z":5, "G":57, "A":71, "S":87, "P":97, "V":99, "T":101, 
             "C":103, "L":113, "N":114, "D":115, "Q":128, "E":129, "M":131, 
             "H":137, "F":147, "R":156, "Y":163, "W":186, "I":113, "K":128}
peptide = PeptideIdentification(Spectrum, Proteome, MassTable)
print(peptide)  # Output: ZXZXX



# Example 2
# -----------
Spectrum = "-10 8 -15 -8 28 29 8 -17 -8 19 6 -9 7 16 -15 -13 -9 29 21 17 -12 28 22 -7 5 -20 -8 13 15 20 -11 21 -5 8 20 21 16 8 -4 4 -20 25 -18 11 15 20 12 -18 -6 -19 -10 -5 10 7 26 -10 -18 18 -8 8 -7 -4 -13 21 -6 3 23 17 24 -10 -13 7 -8 13 -11 13 -12 5 9 -7 18 3 14 25 -20 8 -19 -15 28 14 -3 26 21 23 -18 2 -15 23 -4 0 8 4 10 15 27 1 1 -8 19 10 15 18 -15 22 -5 3 1 4 6 24 16 19 25 -15 23 10 24 23 13 0 -8 -8 28 23 -20 15 -13 13 25 24 -4 20 30 -6 4 28 26 27 1 -20 5 -19 22 -7 -9 10 3 30 13 9 28 4 -7 -16 18 -15 -5 -9 7 29 24 18 8 27 29 8 11 -18 3 2 18 2 0 16 -4 1 -10 -4 18 27 26 -5 28 1 -17 -6 -2 3 -4 10 13 6 29 8 -20 4 16 18 29 9 -8 -17 21 21 -1 2 27 -7 15 27 20 -19 -13 16 12 -19 -20 -11 18 17 19 30 4 -4 9 -4 -9 -13 -6 6 28 2 -11 -2 30 -3 15 19 3 14 21 -14 -1 -17 9 0 7 -3 2 -10 -20 -3 -17 -10 14 8 1 28 24 -6 -15 6 12 27 13 25 26 11 1 -12 16 -14 11 5 -17 21 28 -4 14 13 -11 2 24 28 -7 1 -2 2 14 1 7 17 7 26 5 -18 20 23 6 22 29 13 1 -18 27 21 -20 -6 -7 -12 -6 7 12 17 15 30 11 19 9 7 -10 23 -8 26 15 26 19 -18 -9 12 -5 -10 -4 23 29 11 -7 20 -5 -14 -16 -8 -16 10 9 -13 -4 -13 0 -7 1 4 22 24 4 24 9 15 23 24 17 -10 27 14 13 -3 -12 17 13 1 18 -17 -6 1 15 30 -14 13 27 -18 3 5 20 13 17 26 11 -8 8 25 -9 1 17 -18 -15 -12 -11 16 9 20 19 -2 15 17 -14 2 -4 -14 -17 -19 -12 8 -14 22 -13 1 30 -20 5 28 9 -4 -18 -16 -20 8 24 3 -2 21 27 22 -17 -15 22 20 0 0 8 2 9 -5 -2 -9 16 -8 -17 8 23 8 -17 6 -4 27 -9 -9 -5 -16 -6 -3 -4 16 20 14 30 5 -14 3 23 -3 -14 29 -4 -15 20 -10 13 19 -15 3 5 23 10 -10 -19 8 -17 17 1 -18 4 30 20 -4 -7 14 28 6 27 -11 12 16 4 25 -19 7 -1 13 19 -1 18 -9 15 2 -20 0 15 -1 8 2 11 17 7 -3 -2 -12 -17 -2 -15 22 23 17 -13 21 12 30 7 7 -12 -19 25 23 -1 1 24 20 -16 -13 23 19 14 18 7 -9 10 1 23 27 3 19 12 -4 5 -2 -7 -1 23 -5 -2 25 -20 -16 22 4 -17 -1 14 26 -16 29 -15 3 -5 19 11 1 30 13 -19 4 -15 2 -14 19 -11 -11 1 -16 3 -5 -5 0 23 -3 -10 -3 20 11 24 17 26 11 8 14 -10 -4 -13 -16 25 22 5 22 18 -10 -10 16 0 22 -19 12 26 15 -12 2 5 8 8 -5 -16 9 -19 8 10 -16 -9 -16 -10 -3 -3 -15 1 25 19 -10 21 2 13 9 30 9 -10 25 -16 12 19 6 -12 -9 15 12 -12 0 3 2 15 18 -12 2 7 20 -10 19 14 8 16 0 4 18 0 -5 23 14 -17 6 15 -15 -5 25 25 14 -2 16 -19 -3 13 16 -7 24 10 22 -18 25 2 28 -16 -16 21 4 -2 -16 -20 29 -3 8 9 30 11 6 13 19 26 -1 -7 25 -11 -9 7 16 22 -7 -7 21 0 23 24 -14 26 3 -5 -6 -18 5 16 14 13 19 6 23 1 -1 -12 -18 -5 21 19 26 -5 -8 -19 -8 24 15 11 18 -15 -18 -14 -7 22 -5 18 -4 15 19 4 28 2 -19 1 10 27 -10 -19 27 -10 0 -17 -5 -3 22 -5 30 -1 2 26 -13 -7 -12 -4 19 -8 19 25 15 -1 9 2 -2 27 27 -19 19 14 -10 22 -14 -20 22 23 28 11 -6 -12 26 30 1 27 26 -17 26 -2 8 -20 25 27 -13 14 24 0 -7 1 -5 12 -14 12 21 28 -1 12 19 2 27 -1 8 26 -7 10 30 -6 -18 4 -14 8 27 29 10 -3 9 6 -1 -8 -10 -8 -9 -16 5 26 -13 -16 24 19 -20 -11 4 -12 20 25 4 30 -9 25 -13 18 -6 11 -3 18 29 24 -19 28 18 -14 10 27 27 2 26 1 28 -11 28 24 5 24 6 12 -13 21 -14 17 -6 27 -14 22 20 -16 -12 10 24 25 -10 25 3 -1 -6 -8 -5 -10 -6 -2 -9 29 1 23 12 19 26 -5 -8 -12 -3 15 21 30 -14 -18 11 -2 -13 -20 -19 5 28 -17 8 11 26 -12 7 9 -11 25 11 -11 -18 10 1 -18 12 25 12 -19 -8 23 27 7 23 -10 -12 -15 -6 -9 19 15 14 -11 12 17 -20 9 27 24 -9 -6 -10 15 25 -11 16 -8 26 19 12 13 -4 -6 -6 20 -4 8 1 3 15 23 -10 29 13 -12 -16 -9 14 23 25 16 8 28 25 -20 -8 23 25 26 9 -4 -8 25 11 -12 -11 0 11 3 -5 -18 6 19 -18 18 -6 26 -3 15 19 4 2 -6 3 28 8 22 -2 19 -9 4 8 22 -16 4 -20 -14 -4 28 15 4 23 -15 15 18 -10 -7 19 -18 16 -13 24 -11 -19 8 28 24 -14 -15 11 15 -8 4 -13 3 -4 23 -17 6 6 -7 -8 -15 8 18 -17 13 6 14 17 26 -18 -2 -4 -19 19 9 -9 -2 -1 -9 24 30 18 8 6 20 -10 -8 20 -12 -13 -19 2 -20 -15 11 16 -11 18 11 -20 20 19 -5 22 -2 14 18 10 26 5 26 9 18 30 -8 -5 29 -13 29 19 28 -14 -8 -11 3 -19 17 9 25 7 -18 14 29 -13 4 2 5 -15 -19 -6 15 16 18 18 -3 9 26 9 2 4 -13 6 -9 21 12 14 17 11 25 -11 -18 17 -11 3 16 13 -19 0 -7 7 10 -10 -3 -10 -8 3 -11 -16 -1 18 2 -17 26 16 -10 -18 30 9 24 -14 19 6 -16 -17 -6 10 13 -5 8 16 -5 -17 -4 26 14 -12 4 -4 18 14 7 -7 -15 8 -1 -10 -17 17 -11 -9 -15 10 1 29 -3 20 5 -7 -10 1 30 23 23 -13 5 -7 3 7 2 -8 -6 20 -11 21 0 -17 30 28 17 25 30 14 -5 16 -19 23 18 12 24 8 7 5 3 12 1 -18 -4 -16 27 29 14 11 9 -15 4 -18 -10 22 -18 6 7 -19 -17 -18 26 14 9 20 -6 -6 -6 -6 11 29 -5 19 22 -6 -3 29 1 23 7 -9 16 27 -11 17 -12 -2 -11 -14 15 20 -5 -1 17 -3 27 30 25 -11 23 13 -15 -16 24 -15 6 -13 -18 -9 -10 -8 5 -1 -16 19 -17 30 -20 -12 5 27 28 -1 21 -11 1 -18 19 -2 7 -12 26 24 5 20 20 -10 21 -13 24 19 1 3 27 -8 20 10 -3 4 -8 15 2 1 -16 11 29 8 -2 -20 10 -9 23 -9 30 -11 4 26 17 -2 26 28 20 -16 -9 6 28 11 -15 -9 -3 13 25 13 9 18 9 9 -7 -18 3 26 17 -11 -16 18 -7 3 -1 7 -18 4 0 2 1 9 10 21 4 -16 21 27 -20 -12 28 8 22 15 -4 -7 20 23 23 9 -12 -15 -12 -4 26 4 -5 -15 14 30 13 16 30 12 6 -12 -1 5 19 -13 21 -6 21 30 2 2 30 6 19 -8 25 -16 21 11 -12 14 6 23 20 -19 -15 1 15 16 7 -5 9 28 6 -10 -16 28 -13 8 24 -15 -6 8 -7 -5 -16 -10 24 2 -6 2 -17 19 0 10"
Proteome = "EIFWAPNSFNRWCRFPRKIEICVYQPAGFYHMEMRNVKTDKMHAVDGVITWENVFKAALLYHAAWGGNSGGDEMMACAMFMRWMQGIEAQLHNYGNDMTWEMMAILMMLDRGQHWAEWTGGAWRTEQCNNQEIAENEITQVNSFGTVPACRRQDRVQDYFCQNGRSQPIKAISDCHKMIIDISEMIGWSNLPQLSCLFYQNTWPCLAQVCAHIAQRAQQYTQVHHWAQPDLDAFEYGISDDRTSWQFCCANLFSWVMFSSQFCVEKYKVRDRRILQGLACCTETVQKYQERSCCMHKDEWEVYNDEYIPCMDMFNNKVTMAACRTIVEASEPVLKPPPNKMIEMMDNEHVYVAIVDGFWKLEAARSCFSTRNEKKVAREQRFQECDEWYKRNTGIQEWYIEDPANYPWMAAQLFEKWHCGTQTDNKAQIHMPPFLKVMNIKEDTKWPIEFCCRIVNYCDRPIELEIKFPGTWMLCQMLYSSGWGGHEFEPYLRFCLFYRYGKAIRMHKGDTPTQSLPAFVCFRNGDDLNEGVRPWQRCANIRNEWWPEKPPQEDPSDWGNDACCVCNWCDGDRKVFSESDCCNNDANPNWSEQKIIILRTMQLCDHIRMVYYKAHAWTGMSGHFIYKLFKPWCFFGTFVFAQDVMAPDYYWYHCLGHATKNVVGFFHMVQQSSEIQNCHACASHKRKFHQVEKMSWMCQCPLRCFSNVMPLNMAYYIKTYNVNGTDFYAGITVRSSAWHEQFTSHCSWVFCLQWCMGDPVLFYPSTIGKYMHVRSWTRYREFVVFTGAEMNRGYNYNWAYIVCIIMKNCQAQMCFMKTLHRQQCVQHMPWVFHEIFIADEVGKWYLTTYAAHWDCTRTEAPKTYIIPQQNPIDSQRKPTNCQIGWEMGIRNPLCWFHTSWPHSPTFSEATMIYEFRYNYGRPIDWDWKSQRNTTKRNQMDNGWPWYVISMCMWFYTGQEAVMCTQCAGFQFDRERMIHAGPMLIPIQCKIFNSVPPYNWETYRPGFNTFWRYMGECRMKIGKQTVSLNRVYHWEWCFFGKQIDFIICTNPFTKRDCRLDHNTHSNKMYRNVSNAHTTCFKYKYFYHKNKMCFWLNVEPAKCVAWAEFICWTIVMCFFRVQKEQFPQECWMMIIVTSNLDAFPDNSRQKRIGETDPQQHMGKWHRHARDKGCMPREQTWQDKFCTYDPSAEDKYQDCYLISAFTWQYWMEQCCYKPRHDPWLNKLQPQFLMSIDSPNPCFSCINHPDECKVGYKHVSRGRPMEQHCKMQPIQNAFLPRFALQHNSQHSWWWIANWQDKMSARWPRLLRTHCVIVNAWQACNYNDRHCYHFLTFVMYSGFDQNCELRACVVIQCNRTGDYHNQKVYLMIMANHIMLVPRARWFTPQRCINMIECKKNSEANDQYQRLIEFTPWGVEGHEIPNERILQTRHLLLYCMVYQQVQYAKGHEAAIWWFWDTGRIAYCLCLVLRCYWEPESQTTLLTLGAGNNDKMIHNMTFVETYAMPSVASILEKVQMNPDMIMVNEWHKAVWVWDISTLACFHEHHIHPMSICKLIWMSENRFKNFHCPTYGYITRCAYMAQRDYENKHAIWVGTTKFFKTFGFLGYHCPQRIVIEHLIYRPCSTLWCPNWDKKRCQRRGRQQPYANQLSCIAGNGESRIFWTWPQRQIYLYGFHDNHNQMIYLPTSILHSKSNGYAPQFRQYRLSMETMMMTGHQTKPHDWHNKHIHCCYEHDIFLKFVKEVEMYCFNGNNARREIDFQIKYYKHSAFLFQVNMHLYWFKQCGNMSAGNSFRSYKHESFLLYNRDYMKLLLCCDQCDVWCPYNCDRIITCKAGIIIQMATEYPDSQFAVFERGVDDYRVLGDFDAYDAYNTNHVVWQDAYNMWFYNLGPYAMGRQACKPCRGFHGYVMWWEWMALAGGGSIHSYRSERGCASFFWQPNTCGFYGLAMENYREFLCQMDMVVRNKVSDYNTYMMYDWCKPFANLSMFNCKQPMQKTYHLTFGKLEVDVANYMVLFKLHFGCISCDDERLWDMMGRRQVDNMCKEPRPFGFIYYQQTEFWEGLADVTAQYKLNIGKFTMHQELPVPKDNLMNWECSRWFEECHFTQNLFKLRVHYGFAGSALTKVDNCHYSEIQPFTFGKCQCWSDVLKTLQTMYGLNHMYEDEWEYHILETTFCYGNQNMNPFMAHAWFPYSSWQPNHYRPKWLRHLTNMRSGQDMVTRVVCNCANVLDGIHSRIPFFIPSMKSRYFVTSTMPECPFPQMNSWAQIASTYQYKTPWSGLVQGFENDPFAHLYVWFARNWDMYRQSGRCAMYFVKSMSHIIIMCAVHGDLHPQCVAYNHYKTPNAEDCNRDNPDSHYFNTANMDEWPTSHTTFNMYHLIDMCSDPVANCHMVCVQSPCLNEDTGLRTKTEYGTMAEHTSDQKTYQCCGMENYSTVRYNSVEHFWMTWNRNHKVWEMCNIFYSHTFWYVPWLIHTQYPGRAGSTVKQKWTNIGLAACIVYESYFVEGKVELIRQRFHNPTIRTKNPAFTPCIGFTLYTYMEEEVKMHLGVRNNYEYNMSEPDCAQKLREFQEVEIHRMAKGLIVTNAHHCFSMAFYQHSYGKIHVSKSKGHQCCQPGYGPSYEGPTVPKSRKLNQKGAYGHPPDQNEYHNRPENGYCECNVCWEQNTVDTPSPTDGLHVHHVCVPGWVNSGGTEDEWHWQTYMVMTDHFKLYQAQEYPSIWQRQSCQWFCECDGHAPEKMADMYAWGPDHNFGSEMIRPRRCLPCRKGNTAQEVVGMMQENESWANWGIWPMQSQYRIINMTSEMIDSTKYQQGFMFKGQHCSFQKWSQTESEFKSIYGHETRLAVGFLWWENWIWLREYIRPAGMAGCEDEYFTVYPWEPVASRYHWDWEPKTMDCQSKYPIRWDKDDQKWMFNYLIWPQVNDGPMAQSQSHRREQYKYHLYDTNGCQWEFQTPHKPHEYSIQQWKDWYFYNRMKGISCNFDNVSTLPFMNHQPVFCEHFEGRQAHWCQCIFNCYELGEFDYMQWVMVKYAKIKHNPTHTQEVNVWNCWCTPQCICQWMDFPCWWRCFMTCWETMYPGRMSHRMPTLQFLYIEDHVVLMNYSEGKDHWGCKVHPGDRPHVIIDDNPHDMWVWVKGDMIAYIKMVQNQRVCMALYQLWALKQWYYMKWMNDETLSFWEWCGAPGWFIGGADRQRKIFWCWVIGDHPNPGFSSGAVCYDLIRLYDMITPFDYTAESRQLNKYNSRIEGITFVGQRIQLQMIKMGVPIQQPGTWPMAFDMFWDPIHMQWSWNQKFGFGEIYICHFVGGFAAYMIQNWQIQYVEVLASPRYHQDGIKVHHFLESDQSCTMNNPKIHLDFALGEHFPWTRIHAPSCHFYRNQCLWLVIQHQISQHCWMKCQWAFCLKQTKNILHKESYSVCSAEIAMWTIRGRYLAVRASELIAIIGHELPIQGTYRITANAKDEQCTLRGLKVYYSKKRGEKDVDNCIYEIANGIWTQLNIFVRADWSPMAEWGSLSFYMFPNVDDLYRCHDYVMSYMYYFPQRIYWATCKDWEAKLPRVAENVCRDHMHMGHKRWFHGTEQRRHVLQDQICQKTDRYERWRNDQDQNQMQSGPCYCSSHYFCNEWGQTERFIKYEPKNQSFPRHNYLKRHMKPEEANARSESRYCEMSGQFLWLGLENSKHLVRSERWSRDKYYQKFIARIAAHDKILMINYWRATKGAMSHKRTGPQDGVCIWVCWNSKVRMQGYDCTNNACCWCPSCAAMIKENILNYTHAVNTRVWQFNYITRDYQDCWKERNIRQQPWHRFHSCLQLCCHDFIIDKTTLVDFGSLVCARLHPQCNVVYYKAQINRRFMTDCIDVTIVVWPKQEKNAFGRHHQALFPRGHYSGADCVIFDHTKILRIDLSKPHESHWWKYHLPKIAFAGHSVDDPIFIIYMCNHRVPSWAMAWWVKHLTSSMGVFRCISWCNQHLTRTIRKLMYSFGGAPSLPFQMFGINSLEVPNCVTPFIDYRSRVMSDNHSPLKVDVNIMKWQTEYFKGLWCEMTHMERIDYKLCRQFETWCCDMQHSKMPGQMPFIFVFDEMLIFQTPRPDTYENYYRIIGEMHRNYACWPIAFPHDDLPHTVVLKFANNCFRWHDQVINYFGLFVKVYACLPLWVLLLHCEWMMHRVSRMNVCGPILPCEKADILQVKSHEMFPGANGECPYWDPGSCHDAVACCFAKGGEHCGSEMMNKPANTKCHPSANSMWQRHASPRQIIAANHCIALGKERGYHPLQKKHAWRGLEPQWHKSTYKQRCCPMICADIMVDMTWDAHQWGSNFRDEPIYCWAINTHGLHHEIFCYAILHDHCGEHTMVHWVMFNIPPMQHWKVRHEHPGCRGVMDEPSLVFFSSDGSNWIKVTSEALFVDIYPQQPKGDAWEVCKQEPPEISELGSEHSANASCNNFIYSGTKTFRMDNQGECSNKAFALCFIMCIWCTYTEKYMAPRNAVLSKADCSPIWMEEILFVFMPAYYPGDIFPWMAEKDCIYAMTDNNHDDPVRCECKSEAQEQIMCMPSHFKCRSYVVKRHMSGAPETIGFAWKPCNNGTLMVLKATLAFGYSCGCYLITRSLAQHYMNIYEYIIAMDWCRGMGNEIHPVIDDSCKTRDVTIHCHKDVMWPGEMEENEGCRWVLCVIKQMQVMIFWDINTYHSVNSARIKWRGKFTCGATGHDGWKCCACDEYRIMEFLDLDWQEGGRQMFACVPANSMPFCIFLVVVGYTICHSIICDTNYGYWECPSSECIRVKKSTLYILASHISYKRTTTDLSHYFMMTYYYEFCEVLTMVQTMNSKVVRDYHAEHVQPRVNAVYGPNSYETIVGHERFLTWTCAKWWSSVWNLEMDRRFMLQLFRCEISMQMLFKVQPYYMARACNYGVGYMWEGTATVACLGKVCNVAVFKFEPFCKMLDRPTIQPDIIYNATIYIDDLDSDKICYMSICVYFYNLVSYYSAHEDFAVCIQHITWELDVFGQNDYMELVIGHI"
MassTable = {"G":57, "A":71, "S":87, "P":97, "V":99, "T":101, "C":103, 
             "I":113, "L":113, "N":114, "D":115, "K":128, "Q":128, "E":129, 
             "M":131, "H":137, "F":147, "R":156, "Y":163, "W":186}
peptide = PeptideIdentification(Spectrum, Proteome, MassTable)
print(peptide)  # Output: KLEAARSCFSTRNE



# Example 3
# -----------
with open("dataset_30270_2.txt", "r") as file:
    Spectrum = file.readline().strip()
    Proteome = file.readline().strip()
MassTable = {"G":57, "A":71, "S":87, "P":97, "V":99, "T":101, "C":103, 
             "I":113, "L":113, "N":114, "D":115, "K":128, "Q":128, "E":129, 
             "M":131, "H":137, "F":147, "R":156, "Y":163, "W":186}
peptide = PeptideIdentification(Spectrum, Proteome, MassTable)
print(peptide) # Output: SSISNYSCEVIVICS



# -----------------------------------------------
# Find Peptide-Spectrum Matches (PSMs)
# -----------------------------------------------

def PeptideIdentification(spectrum, proteome, massTable):
    '''Identify all peptide-spectrum matches scoring above a threshold for a set of spectra and a proteome.  
    Input: A set of space-delimited spectral vectors SpectralVectors, an amino acid string Proteome, and an integer threshold.
    Output: The set PSMthreshold(Proteome, SpectralVectors).'''
    
    length = len(proteome)
    max_value = -1000000000
    max_seq = ''
    for i in range(length):
        for j in range(i+1, length):
            if Weight(proteome[i:j+1], massTable) == len(spectrum):
                if PeptideScore(proteome[i:j+1], spectrum, massTable) > max_value:
                    max_value = PeptideScore(proteome[i:j+1], spectrum, massTable)
                    max_seq = proteome[i:j+1]
            if Weight(proteome[i:j+1], massTable) > len(spectrum):
                break

    return max_seq
            
        
def PeptideIdentificationSubs(spectrum, subs, massTable):
    max_value = -1000000000
    max_seq = ''

    for sub in subs:
        if PeptideScore(sub, spectrum, massTable) > max_value:
            max_value = PeptideScore(sub, spectrum, massTable)
            max_seq = sub
    return max_seq, max_value
    
def Weight(seq, massTable):
    total = 0
    for i in seq:
        total += massTable[i]
    return total

def PeptideScore(seq, spectrum, massTable):
    value = 0
    index = -1
    for i in seq:
        index += massTable[i]
        value += int(spectrum[index])
    return value
        
def PSMSearch(spectral_vectors, proteome, threshold, massTable):
    PSMSet = []
    length = len(proteome)
    vector_length = []
    for i in spectral_vectors:
        vector_length.append(len(i.strip().split()))

    sub_weight = dict((i,[]) for i in vector_length)
    for i in range(length):
        for j in range(i+1, length):
            if Weight(proteome[i:j+1], massTable) in vector_length:
                sub_weight[Weight(proteome[i:j+1], massTable)].append(proteome[i:j+1])
            if Weight(proteome[i:j+1], massTable) > max(vector_length):
                break
    
    for vector in spectral_vectors:
        vector = vector.strip().split()
        if len(vector) in sub_weight:
            seq, value = PeptideIdentificationSubs(vector, sub_weight[len(vector)], massTable)
            if value >= threshold:
                PSMSet.append(seq)
    return PSMSet



# Example 1
# ------------
SpectralVectors = ["-1 5 -4 5 3 -1 -4 5 -1 0 0 4 -1 0 1 4 4 4",
                   "-4 2 -2 -4 4 -5 -1 4 -1 2 5 -3 -1 3 2 -3"]
Proteome = "XXXZXZXXZXZXXXZXXZX"
threshold = 5
massTable = {"X":4, "Z":5, "G":57, "A":71, "S":87, "P":97, "V":99, "T":101, 
             "C":103, "L":113, "N":114, "D":115, "Q":128, "E":129, "M":131, 
             "H":137, "F":147, "R":156, "Y":163, "W":186, "I":113, "K":128}
psms = PSMSearch(SpectralVectors, Proteome, threshold, massTable)
for peptide in psms:
    print(peptide)
    
# Output: XZXZ



# Example 2
# ------------
with open("dataset_30270_7.txt") as file:
    lines = file.read().strip().split("\n")
SpectralVectors = lines[:-2]
Proteome = lines[-2]
threshold = int(lines[-1])

MassTable = {"X":4, "Z":5, "G":57, "A":71, "S":87, "P":97, "V":99, "T":101, "C":103, 
             "I":113, "L":113, "N":114, "D":115, "K":128, "Q":128, "E":129, 
             "M":131, "H":137, "F":147, "R":156, "Y":163, "W":186}

psms = PSMSearch(SpectralVectors, Proteome, threshold, MassTable)
for peptide in psms:
    print(peptide)

# Output: DPFCNCVQHAGWISY RSKKVWDRSWLVV HPPISQIQQIPKCD ILRPTILKHLQCTGI KDIHVMRGCQNEW
    # WIHYPAEWCAQEK GAVARYANVNPNVADK PWPRHEMNMWWGS GSQYCMDYLCLEV VSYPNLKSHYNHAEP
    # FFRFELDTGLNLAI TDGGCPPNYYTKHNA SSYCGRACFRQQIHQ NVFPEASCDYAYVPS PCRFMIVVYLKHL
    # RGSNYYPTYHKAY KSLMLCDFHSFLED DDTQIRIPAELSPQ YIKDDQDQAVYMCM LQKVTELRSYFMM
    # MEAVMELAVGQAMY TEALAVTKWNLIKC TFPNCENSFHGVNF



# -----------------------------------------------
# Size of Spectral Dictionary
# -----------------------------------------------

massTable = {}
with open("amino_acid_integer_mass_table.txt", "r") as file:
    for line in file:
        aa, mass = line.strip().split()
        massTable[int(mass)] = aa

def SpectralDictionarySize(spectrum, threshold, max_score):
    '''Computes the number of peptides whose score against the spectral vector
    is at least `threshold` and at most `max_score`.
    
    Input: A spectral vector Spectrum', an integer threshold, and an integer max_score.
    Output: The size of the dictionary Dictionarythreshold(Spectrum').
    '''
    n = len(spectrum)
    masses = list(massTable.keys())

    # Initialize DP table: dp[i][s] = number of peptides ending at position i with score s
    dp = [ [0]*(max_score+1) for _ in range(n+1) ]
    dp[0][0] = 1  # base case: empty peptide

    for i in range(1, n+1):
        for mass in masses:
            if i - mass >= 0:
                for s in range(max_score+1):
                    prev_score = s - spectrum[i-1]
                    if 0 <= prev_score <= max_score:
                        dp[i][s] += dp[i - mass][prev_score]

    # Count peptides with score >= threshold and <= max_score
    total = sum(dp[n][threshold:max_score+1])
    return total



# Example 1
# -----------
spectrum = "4 -3 -2 3 3 -4 5 -3 -1 -1 3 4 1 3"
spectrum = [4, -3, -2, 3, 3, -4, 5, -3, -1, -1, 3, 4, 1, 3]
threshold, max_score = 1, 8

print(SpectralDictionarySize(spectrum, threshold, max_score)) # Output: 0



# Example 2
# -----------
with open("dataset_30266_3.txt", "r") as file:
    lines = file.read().strip().splitlines()
spectrum = list(map(int, lines[0].split()))
threshold = int(lines[1])
max_score = int(lines[2])

print(SpectralDictionarySize(spectrum, threshold, max_score)) # Output: 


# -----------------------------------------------
# Probability of Spectral Dictionary
# -----------------------------------------------

def ProbabilitySpectralDictionary(SpectrumStr, Threshold, MaxScore):
    # Parse input
    Spectrum = list(map(int, SpectrumStr.strip().split()))
    AminoAcids = [57, 71, 87, 97, 99, 101, 103, 113, 114, 115,
                  128, 129, 131, 137, 147, 156, 163, 186]

    m = len(Spectrum)
    
    # Initialize DP table: dp[i][t] = probability of mass i with score t
    dp = [ [0.0]*(MaxScore+1) for _ in range(m+1) ]
    dp[0][0] = 1.0  # base case: mass 0, score 0

    for i in range(m+1):
        for t in range(MaxScore+1):
            if dp[i][t] > 0:
                for aa in AminoAcids:
                    next_i = i + 1  # Each amino acid adds 1 position in the spectrum
                    if next_i <= m:
                        s_i = Spectrum[i]  # Score contribution of this position
                        next_t = t + s_i
                        if 0 <= next_t <= MaxScore:
                            dp[next_i][next_t] += dp[i][t] / 20.0  # divide by 20 for uniform probability

    # Sum probabilities for sequences that reach mass m with score >= Threshold
    final_prob = sum(dp[m][Threshold:])
    return final_prob



# Example 1
# ------------
spectrum_input = "4 -3 -2 3 3 -4 5 -3 -1 -1 3 4 1 3"
threshold, max_score = 1, 8

result = ProbabilitySpectralDictionary(spectrum_input, threshold)
print(f"{result:.3f}")




# -----------------------------------------------
# Spectral Alignment
# -----------------------------------------------

def SpectralAlignment(peptide, spectrum, k):
    n = len(peptide)
    # DP table: dp[i][j] = max score for first i amino acids using j modifications
    dp = [[-float('inf')] * (k + 1) for _ in range(n + 1)]
    back = [[None] * (k + 1) for _ in range(n + 1)]

    dp[0][0] = 0

    for i in range(1, n + 1):
        for j in range(k + 1):
            # Case 1: no modification
            score_no_mod = dp[i-1][j] + spectrum[i-1]
            if score_no_mod > dp[i][j]:
                dp[i][j] = score_no_mod
                back[i][j] = (i-1, j, 0)

            # Case 2: apply a modification if j > 0
            if j > 0:
                # Try reasonable delta modifications
                for delta in range(-10, 11):  # You can adjust range for larger spectra
                    score_mod = dp[i-1][j-1] + spectrum[i-1] + delta
                    if score_mod > dp[i][j]:
                        dp[i][j] = score_mod
                        back[i][j] = (i-1, j-1, delta)

    # Backtracking
    res = []
    i, j = n, k
    while i > 0:
        prev_i, prev_j, delta = back[i][j]
        if delta == 0:
            res.append(peptide[i-1])
        else:
            res.append(f"{peptide[i-1]}({delta:+})")
        i, j = prev_i, prev_j

    return ''.join(res[::-1])




# Example 1
# -----------
peptide = "XXZ"
spectrum = [4, -3, -2, 3, 3, -4, 5, -3, -1, -1, 3, 4, 1, -1]
k = 2

print(SpectralAlignment(peptide, spectrum, k))
# Output: X(+10)X(+10)Z


# Example 2
# -----------
with open("dataset_30269_3.txt", "r") as file:
    peptide = file.readline().strip()
    spectrum = list(map(int, file.readline().strip().split()))
    k = int(file.readline().strip())
            
result = SpectralAlignment(peptide, spectrum, k)
print(result) # Output: M(+10)A(+10)D(+10)SSR

