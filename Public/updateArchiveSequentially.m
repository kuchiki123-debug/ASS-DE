function outArchive = updateArchiveSequentially(archive,pop,fitness)
%% Update the archive sequentially
%writen by zhenyu wang on 20210727
%archive denotes the original archive,and archive consists of NP, pop,
%index and flage.
%pop is the data that will be added into archive.
%fitness denotes the fitness of pop.
    if archive.flag == 0
        outArchive = updateArchive(archive, pop, fitness);
    else
        [px,~] = size(pop);
        for i = 1:px
            archive.pop(archive.index,:) = pop(i,:);
            archive.fitness(archive.index,1) = fitness(i,1);
            archive.index = archive.index + 1;
            if archive.index > archive.NP
                archive.index = 1;
            end
        end
        outArchive = archive;
    end
end

