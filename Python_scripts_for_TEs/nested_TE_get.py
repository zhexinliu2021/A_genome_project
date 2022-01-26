##coding = utf -8

## The script is mainly for mining all the inserted TEs in nested TE sequences !!


Emmer_dir='/home/lab706/jerryliu/Agenome_project/output_RM_Emmer/formated_TE/'
CS_dir='/home/lab706/jerryliu/Agenome_project/CS_format_conversion/formated_TE/'
Tu_dir='/home/lab706/jerryliu/Agenome_project/Tu_format_conversion/formated_TE/'

file_tag=['chr1A','chr2A','chr3A','chr4A','chr5A'\
          ,'chr6A','chr7A','chrUn']

def getnu(string):
    ## ## [ ( id,(xx,xx)), ...]
    return int(string[1][0])

## change the [strat : end ] in the string into '1'.
def change(string, START, END, start, end):
    middle=(end-start +1)*'1'
    left =string[:(start-START)]
    right = string[(end-START + 1):]
    return left + middle + right

def getnumber(string):
    ##123..234
    number=string.split('..')[0]
    return int(number)

def getIDnum(string):
    ##123_Tu
    number = int(string.split('_')[0])
    return number

## transform the number string into positon pairs
def transform(string, START,END):
    list_1=[]
    i=0
    point='off'
    for bp in string:
        if bp == '0' and point =='off':
           start=START+i
           i+= 1
           point = 'ON'


        elif bp =='0' and point == 'ON' and i != (END - START ):
            i+= 1
        elif bp =='0' and point =='ON' and i == (END - START):
            end=START +i
            list_1.append((start, end))
        elif bp =='1' and point == 'ON':
            end=START + i -1
            point= 'off'
            list_1.append((start, end))
            i+=1
        elif bp =='1' and point =='off':
            i+= 1

    return list_1
### [(1, 33), (99, 100)]



spe_list=[]


for tag in file_tag:
    input_file=open(Emmer_dir+'Emmer_'+tag+'_TE.formated','r')

    ##step 1: generate the nested and single TE dict.  sepe_dict={
    ##  /id="1_Tu1"   /post="RLG_famc1.5      /status="fragmented"    C     :  (1,4284)

    sepe_dict={}
    for num, line in enumerate(input_file):
        line=line.strip().split('\t')

        id=str(num)+'_Emmer'
        line[1] = id
        pos_list=[]

        for pair in line[-1].split(','):
            try:
                start=int(pair.split('..')[0])
                end=int(pair.split('..')[1])
            except:
                print('error pair !')
            else:
                pos_list.append(start)
                pos_list.append(end)

        try:
            START=min(pos_list)
            END=max(pos_list)
        except:
            print('empty pos_list')
        else:
            position_array=(START,END)

            name='\t'.join(line[0:5])
            if name not in sepe_dict.keys():
                sepe_dict[name] = position_array
    input_file.close()
    print(len(sepe_dict))

    ##step 2: generate the nested TE dict. nested_dict={
    ##  'chr1A\t21055_Tu\t/post="RLG_famc13\t/status="fragmented"\tC' : { id:(xx,xx), id:(xx,xx),
    ##  id:(xx,xx), ...} , ID:{...}, ...}

    nested_dict={}
    START=0
    END=0
    current_pos=(START, END)
    ID=''
    for id, pair in sorted(sepe_dict.items(), key=getnu):
        '''[ ( id,(xx,xx)), ...]'''
        start=pair[0]
        end=pair[1]
        START=current_pos[0]
        END=current_pos[1]


        if start > END:
            current_pos=(start, end)
            ID=id
            nested_dict[ID] = {}
            if ID not in nested_dict[ID].keys():
                nested_dict[ID][ID] = {}
                nested_dict[ID][ID] = (start, end)

        elif start>= START and end <= END:
            if id not in nested_dict[ID].keys():
                nested_dict[ID][id]={}
                nested_dict[ID][id]=(start, end)


    ##step 3: find the single TE and other nested TEs in the big nested TEs.

    final_total_dict={}
    ## used to carry all the finished TEs


    for nested in nested_dict.keys():
        for query_id, query in nested_dict[nested].items():
            final_pos_list=[]

            ## within each nested TE
            START = query[0]
            END = query[1]

            seq='0'*(END-START+1)


            for sb_id, sb_pair in nested_dict[nested].items():
                start = sb_pair[0]
                end = sb_pair[1]

                if start > START and end< END:
            ## changing the seq within the [start : end] into '1'.
                    seq=change(seq, START, END, start, end)

            ## transform the number string into position pairs
            try:
                pair_list=transform(seq, START, END)
            except:
                print('transform error!')
            else:
                for pair_tu in pair_list:
                    try:
                        start_nu=pair_tu[0]
                        end_nu=pair_tu[1]
                    except:
                        print('error final number')
                    else:
                        final_pos_list.append(str(start_nu)+'..'+str(end_nu))

                final_pos_list=','.join(sorted(final_pos_list, key=getnumber))

                #final_list=query_id.split('\t').append(final_pos_list)

                index=(query_id.split('\t')[1])
                final_list=query_id.split('\t')
                final_list.append(final_pos_list)
                #'\t'.join(final_list)

                final_list='\t'.join(final_list)
                final_total_dict[index] = final_list
                #print(final_list)

    ### outputing the file

    #output_file=open('final_version_TEs2/'+'Emmer_'+tag+'_TE.formated','w')

    #for index in sorted(final_total_dict.keys(), key=getIDnum):
    #    output_file.write(final_total_dict[index]+'\n')
    #    print(final_total_dict[index])
    #output_file.close()






















