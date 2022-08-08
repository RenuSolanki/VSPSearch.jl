using FASTX
using BioSequences

const ids = [
    "NC_014760.1",
    "NZ_CP068731.1",
    "NZ_CP022597.1",
    "NZ_CP022594.1",
    "NZ_CP068732.1"
    ]

const _stop1="TAA"
const _stop2="TAG"
const _startcod = "ATG"
const _sub1 = "TTTATAGCCTTAAAAGGAGAGGATAAATTT"
const _vicinity = 50
#const sub_gs = gs_table(_sub1)
const lensub1 = length(_sub1)

function download_and_save(accession_nums::Vector{String},path::String)
    if !isdir(path)
        println(path, " doesn't exist")
        return
    end
    urlmain = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&"
    urlsuffix = "&rettype=fasta&retmode=text"
    for s in accession_nums
           try
            download(urlmain*"id="*s*urlsuffix, joinpath(path, s*".fasta"))
           catch
            println("Accession number ", s, " could not be downloaded")
           end
      end
end

function download_and_save(accession_nums::String,path::String)
    if !isfile(accession_nums)
        println(accession_nums, " doesn't exist")
        return
    end
    if !isdir(path)
        println(path, " doesn't exist")
        return
    end
    urlmain = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&"
    urlsuffix = "&rettype=fasta&retmode=text"
    open(accession_nums) do f 
        while !eof(f)        
           s = readline(f)
           try
            download(urlmain*"id="*s*urlsuffix, joinpath(path, s*".fasta"))
           catch
            println("Accession number ", s, " could not be downloaded")
           end
        end
    end
end

function substr(s::AbstractString, from::Int, to::Int)
    n = length(s)
    to = min(to, n)
    if (from < 1)||(to<from)||(from>n)
        return ""
    end
    return last(first(s, to), to-from+1)
end

function find_all_occurences(subs::String, s::String, n::Int, vic::Int=_vicinity) # returns positions of first "ATG"
    m = length(subs)
    lst = Int[]
    k = n 
    h = 0 
    for i in 1:(n-m)
        res = findfirst(subs, s)
        if res == nothing
            return lst
        end
        l = res[1]
        j = res[end]
        snew = substr(s, j+1, k)
        res1 = findfirst(_startcod, snew)

        if (res1==nothing)
            return lst
        elseif (res1[1]>vic)
            s=snew
            k = k-j
            h +=j 
            continue
        end

        k = k-j #length of new s 
        h +=j 

        push!(lst, h+res1[1])
        
        s = snew
    end  
    return lst
end


function find_all_occurences(subs::Nothing, s::String, n::Int, vic::Int=_vicinity) # returns positions of first "ATG"
    lst = Int[]
    k = n 
    h = 0 
    for i in 1:(n-3)
        res = findfirst(_startcod, s)

        if (res==nothing)
            return lst
        else
            s= substr(s, res[1]+3, k)
            k=k-(res[1]+2)
            push!(lst, h+res[1])
            h+=(res[1]+2)
        end
    end  
    return lst
end

function find_orf(dna, endpos,n, stops::Vector{String}=[_stop1])
    lastpos = endpos
    found = false
    stopcod = ""
    for i in 1:(Int(ceil((n-endpos)/3))+1)
        str = substr(dna, lastpos, lastpos+2)
        if (str in stops)
            found = true
            stopcod=str
            break
        end
        lastpos += 3
    end
    if (found)&&(lastpos-3>endpos)
        return substr(dna, endpos, lastpos-1)*stopcod
    else
        return ""
    end

end


function _search_marked_orfs(path::String, files::Vector{String}, marker_seq::Union{String, Nothing}, 
    stops::Vector{String}=[_stop1], vic::Int=_vicinity, ncbi_gencode_table::Int=4)

    ref = 1
    vec = Tuple{String, String}[]

    for fl in files
        reader = FASTA.Reader(open(joinpath(path, fl), "r"))
        id=""
        ds=""
        dna = ""
        for record in reader
            id=identifier(record)
            ds = description(record)
            dna = string(sequence(record)) 
            break #considers only first sequence
        end
        close(reader)

        n = length(dna)

        #dna_rev = string(reverse_complement(LongDNA{4}(dna))) #use 2 instead of 4 if there are no ambiguities in the sequence; it probably doesn't matter in our case.

        vecpos = find_all_occurences(marker_seq, dna, n, vic) # positions of first "ATG" after conserved sequence


        if length(vecpos)>0
            num = 1
            for pos in vecpos
                orf = find_orf(dna, pos,n, stops)
                lorf = length(orf)
                if lorf>0
                    last_pos = pos + lorf-1
                    name = id#*" "*ds
                    push!(vec, (name*":"*"orf"*string(num)*"-"*string(pos)*"..."*string(last_pos), orf))
                    num+=1
                end
            end
        end
        println(fl, " search done by thread ", Threads.threadid())

        ref+=1
    end

    w = open(FASTA.Writer, joinpath(@get_scratch!("temp_data_dir"), "outfile"*string(Threads.threadid())*".fasta"))
    for x in vec
        des = x[1]
        dna = LongDNA{4}(x[2])
        rna=convert(LongSequence{RNAAlphabet{4}}, dna)
        pro = translate(rna, code=ncbi_trans_table[ncbi_gencode_table])
        rec = FASTA.Record(des, pro)
        write(w, rec)
    end
    close(w)

end



function search_marked_orfs(path::String; marker_seq::Union{String, Nothing}=nothing, 
    stops::Vector{String}=[_stop1], vic::Int=_vicinity, ncbi_gencode_table::Int=4,
    outfile::String=joinpath(@get_scratch!("search_data"), "outfile.fasta"))

    if !isdir(path)
        println(path, " doesn't exist")
        return
    end
    files = readdir(path)
    len_files = length(files)
    println(len_files, " number of ids to be processed ---")

    num_thr = Threads.nthreads()

    m = Int(ceil(len_files/num_thr))

    files_vec = Vector{Vector{String}}()
    for j in 1:num_thr
        a = 1
        b = m
        a = m+1
        b = 2m
        a = 2m+1
        b = 3m
        if j<num_thr
            a=(j-1)*m+1
            b = j*m
            strs = b <= len_files ? files[a:b] : String[]
            push!(files_vec, strs)
        else
            a = (j-1)*m+1
            b =  len_files
            strs = files[a:b]
            push!(files_vec, strs)
        end
    end

    Threads.@threads for i in 1:num_thr
        _search_marked_orfs(path, files_vec[i], marker_seq, stops, vic, ncbi_gencode_table)
    end
    
    temp_data_dir = @get_scratch!("temp_data_dir")
    temp_data_files = readdir(temp_data_dir)

    try
        w = open(FASTA.Writer, outfile)

        for fl in temp_data_files 
            reader = FASTA.Reader(open(joinpath(temp_data_dir, fl), "r"))
            for rec in reader
                write(w, rec)
                #rec = FASTA.Record(x[1], LongDNA{4}(x[2]))
            end
            close(reader)

        end
        close(w)
        println("Orf data saved successfully at ", outfile)
    catch
        println("Could not write data to file ", outfile)
    end

    rm(temp_data_dir, recursive=true)


end


function run_blastp(;infile::String, outfile::String, database::String, other_arguments::Vector{String}=String[])
    command = `blastp -query $infile -db $database -out $outfile $other_arguments` #can also set -culling_limit 2 https://www.ncbi.nlm.nih.gov/books/NBK279684/
    run(command)
    return
end

function run_cdhit(;infile::String, outfile::String, other_arguments::Vector{String}=String[])
    command = `cd-hit -i $infile -o $outfile $other_arguments`
    run(command)
    return
end


function vsp_az(;accession_nums::Union{String, Vector{String}},
    download_folder::String,
    blastpdb::String,
    marker_seq::Union{String, Nothing}=nothing, 
    stops::Vector{String}=[_stop1], 
    vic::Int=_vicinity, 
    ncbi_gencode_table::Int=4,
    orfoutfile::String=joinpath(@get_scratch!("search_data"), "outfile.fasta"),
    cdhit_outfile::String =joinpath(@get_scratch!("cd_hit_res"), "cdhit_out"),
    cdhit_extra_args::Vector{String} = String[],
    blastp_extra_args::Vector{String} =String[],
    blastp_outfile::String=joinpath(@get_scratch!("blastp_res"), "blastp_out")
    
    )

    println("Downloading.....")
    download_and_save(accession_nums, download_folder)
    println("Searching orfs....")
    search_marked_orfs(download_folder,marker_seq=marker_seq,stops=stops,
    vic=vic,ncbi_gencode_table=ncbi_gencode_table, outfile=orfoutfile)
    println("Running cd-hit....")
    run_cdhit(infile=orfoutfile,outfile=cdhit_outfile, other_arguments=cdhit_extra_args)
    println("Running blastp....")
    run_blastp(infile=cdhit_outfile, outfile=blastp_outfile, database=blastpdb, other_arguments=blastp_extra_args)

    println("blast result is saved at ", blastp_outfile)

    
end







    