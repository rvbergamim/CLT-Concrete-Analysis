function Titulo(fid, titulo, larguraTotal)
    % Calcula o número de ":" a serem adicionados antes e depois do título
    numPontos = (larguraTotal - length(titulo)) / 2;

    % Gera a string de ":" para o início e o final do título
    linhaPontos = repmat(':', 1, floor(numPontos));

    % Monta a linha do título centralizado
    linhaTitulo = [linhaPontos titulo linhaPontos];

    % Exibe o título centralizado
    %disp(linhaTitulo);
    fprintf(fid,'\n');
    fprintf(fid,linhaTitulo);
    fprintf(fid,'\n');
end
