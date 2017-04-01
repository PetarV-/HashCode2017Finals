#include <iostream>

using namespace std;

const int N = 1000;
int n, m, radius, cost_edge, cost_router, budget, start_i, start_j;
char board[N][N];

int main()
{
    scanf("%d %d %d", &n, &m, &radius);
    scanf("%d %d %d", &cost_edge, &cost_router, &budget);
    scanf("%d %d", &start_i, &start_j);

    for(int i = 0; i < n; i++)
	for(int j = 0; j < m; j++)
	    scanf(" %c", &board[i][j]);

    printf("dimensions: %d x %d\n", n, m);
    int total_target = 0;
    for(int i = 0; i < n; i++)
	for(int j = 0; j < m; j++)
	    total_target += board[i][j] == '.';
    printf("total targets: %d\n", total_target);
    int area = (2 * radius + 1) * (2 * radius + 1);
    printf("routers to cover: %d\n", (total_target + area - 1) / area);
    
    return 0;
}
