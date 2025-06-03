function x = gauss_sem_pivoteamento(A_original, B_original)
    clc();
    printf("*** MÉTODO DIRETO: GAUSS (ELIMINAÇÃO GAUSSIANA) SEM PIVOTEAMENTO ***\n\n");

    // Exibe entradas
    printf("Entrada - Matriz A (original):\n");
    disp(A_original);
    printf("Entrada - Vetor B (original):\n");
    disp(B_original);

    A_copia = A_original;
    B_copia = B_original;
    n = length(B_copia);

    // Eliminação
    for k = 1:n-1
        if A_copia(k,k) == 0 then
            error("Pivô nulo encontrado. Método sem pivoteamento falha");
        end
        for i = k+1:n
            m = A_copia(i,k) / A_copia(k,k);
            A_copia(i,k) = 0;
            for j = k+1:n
                A_copia(i,j) = A_copia(i,j) - m * A_copia(k,j);
            end
            B_copia(i) = B_copia(i) - m * B_copia(k);
        end
    end

    // Substituição regressiva
    x = zeros(n,1);
    x(n) = B_copia(n) / A_copia(n,n);
    for k = n-1:-1:1
        soma = 0;
        for j = k+1:n
            soma = soma + A_copia(k,j) * x(j);
        end
        x(k) = (B_copia(k) - soma) / A_copia(k,k);
    end

    // Saída da solução
    printf("\nSolução X do Sistema:\n");
    for i = 1:n
        mprintf(" x(%d) = %.6f\n", i, x(i));
    end

    // Verificação
    printf("\nVerificação dos resultados (AX = B):\n");
    for i = 1:n
        s = 0;
        for j = 1:n
            s = s + A_original(i,j) * x(j);
            if j < n then
                printf(" (%d*%.6f) + ", A_original(i,j), x(j));
            else
                printf(" (%d*%.6f) = ", A_original(i,j), x(j));
                printf(" %.6f\n", s);
            end
        end
    end

    // Erro
    printf("\nErro absoluto (AX - B):\n");
    erro = A_original * x - B_original;
    disp(erro);

    printf("\n ********** ELIMINAÇÃO GAUSSIANA FINALIZADA **********\n");
endfunction

function X = gauss_seidel_guloso(A_original, B_original, epsilon, Nmax)
    clc();
    printf("*** MÉTODO ITERATIVO: GAUSS-SEIDEL (REORDENAÇÃO GULOSA) ***\n");

    n = size(A_original,1);
    X0 = zeros(n,1);
    X = X0;
    
    printf("\n Matriz A original:\n");
    disp(A_original);
    printf("\n Vetor B original:\n");
    disp(B_original);

    // Subfunção de reordenação gulosa
    function [A_greedy, B_greedy, sucesso, ordem] = reordenar_greedy(A, B)
        n = size(A,1);
        usados = zeros(n,1);
        ordem = zeros(n,1);
        sucesso = %T;
        
        for j = 1:n
            maior = -%inf;
            linha_melhor = -1;
            for i = 1:n
                if usados(i) == 0 then
                    if abs(A(i,j)) > maior then
                        maior = abs(A(i,j));
                        linha_melhor = i;
                    end
                end
            end
            if linha_melhor == -1 then
                sucesso = %F;
                A_greedy = A;
                B_greedy = B;
                return;
            end
            ordem(j) = linha_melhor;
            usados(linha_melhor) = 1;
        end
        A_greedy = A(ordem, :);
        B_greedy = B(ordem);
    endfunction

    [A, B, sucesso, ordem_linhas] = reordenar_greedy(A_original, B_original);

    if sucesso then
        printf("\n Reordenação Gulosa aplicada com sucesso.\n");
        printf(" Ordem das linhas escolhida:\n");
        disp(ordem_linhas');
        printf("\n Matriz A após reordenação:\n");
        disp(A);
        printf("\n Vetor B após ordenação:\n");
        disp(B);
    else
        error(" Não foi possível aplicar a reordenação gulosa.");
    end

    // Método de Gauss-Seidel
    for k = 1:Nmax
        for i = 1:n
            S1 = 0; S2 = 0;
            for j = 1:i-1
                S1 = S1 + A(i,j) * X(j);
            end
            for j = i+1:n
                S2 = S2 + A(i,j) * X0(j);
            end
            X(i) = (B(i) - S1 - S2) / A(i,i);
        end
        erro = max(abs(X - X0));
        if erro < epsilon then
            break;
        end
        X0 = X;
    end

    printf("\n Número de iterações: %d\n", k);
    printf("\n Vetor solução aproximada:\n");
    mprintf(" x(%d) = %.6f\n", [(1:n)', X]);

    printf("\n Verificação dos resultados (A*X ≈ B):\n");
    for i = 1:n
        s = 0;
        for j = 1:n
            s = s + A(i,j) * X(j);
            if j < n then
                printf(" (%.1f*%.6f) +", A(i,j), X(j));
            else
                printf(" (%.1f*%.6f) = ", A(i,j), X(j));
                printf("%.6f\n", s);
            end
        end
    end

    printf("\n***** ENCERRAMENTO DO GAUSS-SEIDEL COM MÉTODO GULOSO *****\n");
endfunction

function X = tdma_thomas(a, b, c, d)
    clc();
    printf("*** MÉTODO DIRETO: THOMAS (TDMA) - SISTEMAS TRIDIAGONAIS ***\n");

    // Guarda cópias dos vetores originais para cálculo do erro
    a_original = a;
    b_original = b;
    c_original = c;
    d_original = d;

    printf("\n Vetor a^T:"); disp(a');
    printf("\n Vetor b^T:"); disp(b');
    printf("\n Vetor c^T:"); disp(c');
    printf("\n Vetor d^T:"); disp(d');

    n = length(b);
    if b(1) == 0 then
        error("Pivô nulo encontrado na primeira linha. Método falha.");
    end

    // Etapa de eliminação direta
    c(1) = c(1) / b(1);
    d(1) = d(1) / b(1);

    for i = 2:n-1
        temp = b(i) - a(i) * c(i-1);
        if temp == 0 then
            error("Pivô nulo encontrado na linha " + string(i) + ". Método falha.");
        end
        c(i) = c(i) / temp;
        d(i) = (d(i) - a(i) * d(i-1)) / temp;
    end

    temp = b(n) - a(n) * c(n-1);
    if temp == 0 then
        error("Pivô nulo encontrado na última linha. Método falha.");
    end
    d(n) = (d(n) - a(n) * d(n-1)) / temp;

    // Substituição regressiva
    X = zeros(n,1);
    X(n) = d(n);
    for i = n-1:-1:1
        X(i) = d(i) - c(i) * X(i+1);
    end

    printf("\n Solução X do sistema:\n");
    mprintf("  x(%d) = %.6f\n", [(1:n)', X]);

    // Construção da matriz A para verificação
    A = diag(b_original);
    for i = 2:n
        A(i, i-1) = a_original(i);
    end
    for i = 1:n-1
        A(i, i+1) = c_original(i);
    end

    printf("\n Erro absoluto (A*X - d):\n");
    erro = A * X - d_original;
    disp(erro);

    printf("\n ********** TDMA FINALIZADO **********\n");
endfunction

function X = lu_crout(A, B)
    clc();
    printf("\n***** MÉTODO DIRETO: FATORAÇÃO LU por CROUT *****\n\n")

    // Armazena cópias originais
    A_original = A;
    B_original = B;

    printf("Entrada - Matriz A (original):\n")
    disp(A_original)

    printf("\nEntrada - Vetor B (original):\n")
    disp(B_original)

    n = length(B);
    L = zeros(n, n);
    U = eye(n, n); // U com diagonal principal 1 (Crout)

    // Fatoração LU
    for j = 1:n
        for i = j:n
            soma = 0;
            for k = 1:j-1
                soma = soma + L(i,k) * U(k,j);
            end
            L(i,j) = A(i,j) - soma;
        end

        if L(j,j) == 0 then
            error("Pivô nulo encontrado. Método sem pivoteamento falha.");
        end

        for i = j+1:n
            soma = 0;
            for k = 1:j-1
                soma = soma + L(j,k) * U(k,i);
            end
            U(j,i) = (A(j,i) - soma) / L(j,j);
        end
    end

    // Resolvendo LY = B (substituição progressiva)
    Y = zeros(n,1);
    for i = 1:n
        soma = 0;
        for j = 1:i-1
            soma = soma + L(i,j) * Y(j);
        end
        Y(i) = (B(i) - soma) / L(i,i);
    end

    printf("\nSolução Y (LY = B):\n")
    disp(Y)

    // Resolvendo UX = Y (substituição regressiva)
    X = zeros(n,1);
    for i = n:-1:1
        soma = 0;
        for j = i+1:n
            soma = soma + U(i,j) * X(j);
        end
        X(i) = (Y(i) - soma) / U(i,i); // U(i,i) = 1 para Crout, mas mantido para generalidade
    end

    printf("\nSolução X (UX = Y):\n")
    mprintf(" x(%d) = %.6f\n", [(1:n)', X]);

    // Verificação AX ≈ B
    printf("\nVerificação dos resultados (A*X ≈ B):\n")
    for i = 1:n
        s = 0;
        for j = 1:n
            s = s + A_original(i,j) * X(j);
            if j < n then
                printf("(%d*%.6f) + ", A_original(i,j), X(j));
            else
                printf("(%d*%.6f) = ", A_original(i,j), X(j));
                printf("%.6f\n", s);
            end
        end
    end

    printf("\nErro absoluto (AX - B):\n")
    erro = A_original * X - B_original;
    disp(erro)

    printf("\n ********** FATORAÇÃO LU FINALIZADA **********\n")
endfunction

function [X, k] = gauss_jacobi_guloso(A, B, epsilon, Nmax)
    clc();
    printf("*** MÉTODO ITERATIVO: GAUSS-JACOBI (REORDENAÇÃO GULOSA) ***\n");

    n = size(A,1);
    X0 = zeros(n,1);
    X = X0;
    A_original = A;
    B_original = B;

    printf("\n Matriz A original:\n");
    disp(A_original);
    printf("\n Vetor B original:\n");
    disp(B_original);

    function [A_greedy, B_greedy, sucesso, ordem] = reordenar_greedy(A, B)
        n = size(A,1);
        usados = zeros(n,1);
        ordem = zeros(n,1);
        sucesso = %T;
        
        for j = 1:n
            maior = -%inf;
            linha_melhor = -1;
            for i = 1:n
                if usados(i) == 0 then
                    if abs(A(i,j)) > maior then
                        maior = abs(A(i,j));
                        linha_melhor = i;
                    end
                end
            end
            if linha_melhor == -1 then
                sucesso = %F;
                A_greedy = A;
                B_greedy = B;
                return;
            end
            ordem(j) = linha_melhor;
            usados(linha_melhor) = 1;
        end
        
        A_greedy = A(ordem, :);
        B_greedy = B(ordem);
    endfunction

    [A, B, sucesso, ordem_linhas] = reordenar_greedy(A, B);

    if sucesso then
        printf("\n Reordenação Gulosa aplicada com sucesso.\n");
        printf(" Ordem das linhas escolhida:\n");
        disp(ordem_linhas');
        printf("\n Matriz A após reordenação:\n");
        disp(A);
        printf("\n Vetor B após reordenação:\n");
        disp(B);
    else
        error(" Não foi possível aplicar a reordenação gulosa.");
    end

    // Iterações de Gauss-Jacobi
    for k = 1:Nmax
        for i = 1:n
            S = 0;
            for j = 1:n
                if i <> j then
                    S = S + A(i,j) * X0(j);
                end
            end
            X(i) = (B(i) - S) / A(i,i);
        end
        erro = max(abs(X - X0));
        if erro < epsilon then
            break;
        end
        X0 = X;
    end

    printf("\n Número de iterações: %d\n", k);
    printf("\n Vetor solução aproximada:\n");
    mprintf(" x(%d) = %.6f\n", [(1:n)', X]);

    // Verificação A*X ≈ B
    printf("\n Verificação dos resultados (A*X ≈ B):\n");
    for i = 1:n
        s = 0;
        for j = 1:n
            s = s + A(i,j) * X(j);
            if j < n then
                printf(" (%.1f*%.6f) +", A(i,j), X(j));
            else
                printf(" (%.1f*%.6f) = ", A(i,j), X(j));
                printf("%.6f\n", s);
            end
        end
    end

    printf("\n***** ENCERRAMENTO DO GAUSS-JACOBI COM MÉTODO GULOSO *****\n");
endfunction

// Dados de entrada
//A = [0.1, 0.2, 1.0, 0.3;
//     0.3, 2.0, -0.3, -0.9;
//     4.0, 2.0, -0.3, 0.8;
//     0.6, 3.2, -1.8, 0.4];
//B = [4.0; 7.5; 4.4; 10.0];

// Parâmetros de controle
//epsilon = 1.0e-6;
//Nmax = 100;

// Chamada da função
//[X, iteracoes] = gauss_jacobi_guloso(A, B, epsilon, Nmax);




// Definindo entrada
//A = [3, 2, 4; 1, 1, 2; 4, 3, -2];
//B = [1; 2; 3];

// Chamando a função
//X = lu_crout(A, B);





// Entrada dos dados
//a = [0; -5; -5; -5];
//b = [20; 15; 15; 19];
//c = [-5; -5; -5; 0];
//d = [1100; 100; 100; 100];

// Chamada da função
//X = tdma_thomas(a, b, c, d);




// Dados do sistema
//A = [0.1, 0.2, 1.0, 0.3;
//     0.3, 2.0, -0.3, -0.9;
//     4.0, 2.0, -0.3, 0.8;
//     0.6, 3.2, -1.8, 0.4];
//B = [4.0; 7.5; 4.4; 10.0];

// Parâmetros
//epsilon = 1.0e-6;
//Nmax = 100;

// Chamada da função
//X = gauss_seidel_guloso(A, B, epsilon, Nmax);





//A = [3, 2, 4; 1, 1, 2; 4, 3, -2];
//B = [1; 2; 3];

//x = gauss_sem_pivoteamento(A, B);
