##coding = utf - 8
## The script is mainly for calculating the coverage of those nasted TE!!

## Working directory : /public-supool/home/xlshi/jerryliu/Agenome_project/Emmer
import re


def getnumber(string):
    'chr1A_122'
    tag = int(string.split('_')[1])
    return tag


def getposnum(string):
    '123..4565'
    posnumb = int(string.split('..')[0])
    return posnumb


def getsupername(string):
    regexp = re.compile('["_]')
    supername = regexp.split(string)[1]
    return supername


filename_taglist = ['chr1A', 'chr2A', 'chr3A', 'chr4A', 'chr5A', 'chr6A', 'chr7A', 'chrUn']

for tag in filename_taglist:

    input_file = open('Emmer_' + tag + '_TE.formated', 'r')

    single_TE_dict = {}
    parced_TE_dict = {}
    fam_dict = {}

    '''step 1: generating the three dict '''
    for number, line in enumerate(input_file, start=1):

        line = line.strip().split('\t')
        line[1] = 'ID=chr1A_' + str(number)

        fam_dict['chr1A_' + str(number)] = '\t'.join(line[1:5])
        ''' /post="DTC_famc6.1 \t /status="fragmented"  '''

        try:
            if len(line[-1].strip().split(',')) > 1:

                parced_TE_dict['chr1A_' + str(number)] = []
                a = parced_TE_dict['chr1A_' + str(number)]
                a.extend(line[-1].split(','))

            elif len(line[-1].strip().split(',')) == 1:

                single_TE_dict['chr1A_' + str(number)] = line[-1].strip()
                '''single_TE_dict = { 'chr1A_1' :'123..123',....}'''




        except:
            print(line)

    input_file.close()

    ''' step 2: re-editing each TE in parced_TE_dict '''

    for TE in sorted(parced_TE_dict.keys(), key=getnumber):
        empty_list = []
        ''' each parced TE !!!!'''
        for pair in parced_TE_dict[TE]:
            number1 = int(pair.split('..')[0])
            number2 = int(pair.split('..')[1])

            empty_list.append(number1)
            empty_list.append(number2)

        HEAD = min(empty_list)
        TAIL = max(empty_list)
        pos_list = []

        pos_list.append(str(HEAD) + '..' + str(TAIL))

        # print(pos_list)

        ''' pos_list=['1234..56789'] '''
        ''' ***starting go through all the single_Te_dict and parced_TE_dict'''
        for single in single_TE_dict.keys():
            start = int(single_TE_dict[single].split('..')[0])
            end = int(single_TE_dict[single].split('..')[1])
            '''start  the begaining number of 'single'  '''

            new_list = []

            if start > HEAD and end < TAIL:
                # print( 'got one candidate single TE to be serched of')

                for index, part in enumerate(pos_list):
                    START = int(part.split('..')[0])
                    END = int(part.split('..')[1])

                    if start < START and END >= end >= START:
                        START = (end + 1)
                        END = END

                        new_list.append(str(START) + '..' + str(END))
                        # print('got one ')
                        # print(pos_list)
                    elif start >= START and end <= END:
                        START1 = START
                        END1 = start - 1
                        START2 = end + 1
                        END2 = END

                        new_list.append(str(START1) + '..' + str(END1))
                        new_list.append(str(str(START2) + '..' + str(END2)))
                        # print('got one ')
                        # print(pos_list)

                    elif START <= start <= END and end > END:
                        START = START
                        END = start - 1

                        new_list.append(str(START) + '..' + str(END))




                    else:

                        new_list.append(part)

                        # print('not in the TE,so bad~~~~',pos_list)

                    ''' at this time, the pos_lsit is like this: 
                    [['11..22','33..43'],['123..432'],.....]
                         |                     |
                    [  '123..'321'    ,      '323..532'......]
                    '''

                pos_list = new_list
                ####  each single TE going through has a new pos_list.

        # print('The first round searching: ',pos_list)

        # print('The single_dict searching is finished!!')

        new_list = []

        for parced in parced_TE_dict.keys():

            empty_list = []

            for pair in parced_TE_dict[parced]:
                number1 = int(pair.split('..')[0])
                number2 = int(pair.split('..')[1])

                empty_list.append(number1)
                empty_list.append(number2)

            start = min(empty_list)
            end = max(empty_list)

            new_list = []

            if start > HEAD and end < TAIL:
                # print('The parced TE seatching is:' + str(start) + '..' + str(end))
                for index, part in enumerate(pos_list):
                    START = int(part.split('..')[0])
                    END = int(part.split('..')[1])

                    if start < START and START <= end <= END:
                        START = end + 1
                        END = END
                        # print('conditionA: '+str(START)+'..'+str(END))

                        new_list.append(str(START) + '..' + str(END))

                    elif start >= START and end <= END:
                        START1 = START
                        END1 = start - 1
                        START2 = end + 1
                        END2 = END

                        # pos_list[index] = []
                        new_list.append(str(START1) + '..' + str(END1))
                        new_list.append(str(str(START2) + '..' + str(END2)))
                        # print('conditionB: ' + str(START1) + '..' + str(END1)\
                        # +' and '+str(START2) + '..' + str(END2))


                    elif START <= start <= END and end > END:
                        START = START
                        END = start - 1

                        # pos_list[index] = []
                        new_list.append(str(START) + '..' + str(END))

                        # print('conditionC: '+str(START) + '..' + str(END))
                    # elif start <= START and end >= END:

                    else:

                        new_list.append(part)
                        # print('else :'+part )

                # print('after one parced searching: ',new_list)
                pos_list = new_list

        # print( 'FINAL  : ',pos_list)
        '''step 3 : give the TE a new pos_list:
        pos_list = ['123..345','543..987',...]
        TE = 'chr1A_12' ...
        fam_dict={ 'chr1A_1':'/post="DTC_famc6.1 ^I /status="fragmented"',.....}

        '''
        pos_list = sorted(pos_list, key=getposnum)

        position = ','.join(pos_list)

        fam = fam_dict[TE].split('\t')
        fam.append(position)

        fam.insert(0, tag)
        fam = '\t'.join(fam)

        fam_dict[TE] = fam

    '''  step 4: reforming the fam_dict, adding the positon of single TEs in the dictionary '''

    for single_TE in fam_dict.keys():
        if single_TE in single_TE_dict.keys():
            fam = fam_dict[single_TE]
            fam = fam.split('\t').insert(0, tag)
            fam = fam.insert(-1, single_TE[single_TE])
            fam = '\t'.join(fam)

            fam_dict[single_TE] = fam

    print('the ' + tag + 'is finised')

    ''' step 5: outputing the dict '''
    output_file = open('final_version_TEs/' + 'Emmer_' + tag + '_TE.formated', 'w')

    for TE in sorted(fam_dict.keys(), key=getnumber):
        output_file.write(fam_dict[TE] + '\n')

    output_file.close()

    ''' step 6: calculating the TE coverage '''
    superfam_dict = {}
    '''subfam_dict = { 'DTT':123,'RLG':234,....}'''

    for TE in sorted(fam_dict.keys(), key=getnumber):
        superfam = getsupername(fam_dict[TE].split('\t')[2].strip())

        position_list = fam_dict[TE].split('\t')[-1].split(',')

        if superfam not in superfam_dict.keys():
            superfam_dict[superfam] = 0

        length = 0
        try:
            for pair in position_list:
                start = int(pair.split('..')[0])
                end = int(pair.split('..')[1])

                length += end - start + 1

            superfam_dict[superfam] += length
        except:
            print(fam_dict[TE])

    ### calculating the length of whole genome
    ###  /home/lab706/jerryliu/Agenome_project/output_RM_Emmer

    seq_input = open('/public-supool/home/xlshi/jerryliu/Agenome_project/Emmer/genome/' + tag + '.fasta')

    seq = 0
    for line in seq_input:
        if line[0] != '>':
            seq += len(line.strip())

    coverage = 0
    output_file = open('result', 'w')
    for supername in superfam_dict.keys():
        coverage += round(superfam_dict[supername] / seq, 1)

    output_file.write('The percentage of TE in' + tag + ' is: ' + round(coverage, 1) + '\n')
