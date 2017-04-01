#include <iostream>

using namespace std;

const int N = 1000;
int n, m, radius, cost_edge, cost_router, budget, start_i, start_j;
char board[N][N];

int value(vector<coord> routers)
{
    bool mark[n][m];
    for (int i=0;i<n;i++)
    {
        for (int j=0;j<m;j++)
        {
            mark[i][j] = false;
        }
    }

    int ret = 0;

    for (int i=0;i<routers.size();i++)
    {
        // check if outta bounds
        if (routers[i].i < 0 || routers[i].i > n - 1 || routers[i].j < 0 || routers[i].j > m - 1) return -2;
        // check if wallhacked
        if (mat[routers[i].i][routers[i].j] == '#') return -1;
        // otherwise, expand up
        for (int dx=0;dx>=-radius;dx--)
        {
            int xt = routers[i].i + dx;
            if (xt < 0) break;
            
            // top-left
            for (int dy=0;dy>=-radius;dy--)
            {
                int yt = routers[i].j + dy;
                if (yt < 0) break;
                if (mat[xt][yt] == '#') break; // no need to go on
                if (!mark[xt][yt])
                {
                    mark[xt][yt] = true;
                    if (mat[xt][yt] == '.') ret++;
                }
            }

            // top-right
            for (int dy=0;dy<=radius;dy++)
            {
                int yt = routers[i].j + dy;
                if (yt > m - 1) break;
                if (mat[xt][yt] == '#') break;
                if (!mark[xt][yt])
                {
                    mark[xt][yt] = true;
                    if (mat[xt][yt] == '.') ret++;
                }
            }
        }

        // now expand down
        for (int dx=0;dx<=radius;dx++)
        {
            int xt = routers[i].i + dx;
            if (xt > n - 1) break;
            
            // bottom-left
            for (int dy=0;dy>=-radius;dy--)
            {
                int yt = routers[i].j + dy;
                if (yt < 0) break;
                if (mat[xt][yt] == '#') break; // no need to go on
                if (!mark[xt][yt])
                {
                    mark[xt][yt] = true;
                    if (mat[xt][yt] == '.') ret++;
                }
            }

            // top-right
            for (int dy=0;dy<=radius;dy++)
            {
                int yt = routers[i].j + dy;
                if (yt > m - 1) break;
                if (mat[xt][yt] == '#') break;
                if (!mark[xt][yt])
                {
                    mark[xt][yt] = true;
                    if (mat[xt][yt] == '.') ret++;
                }
            }
        }
    }
    return ret;
}

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
