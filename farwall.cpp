#include <queue>
#include <cstdio>
#include <cassert>

using namespace std;

char test_files[4][50]={"tests/charleston_road.in","tests/lets_go_higher.in","tests/opera.in","tests/rue_de_londres.in"};

const int N = 1000;
int n, m, radius, cost_edge, cost_router, budget, start_i, start_j;
char board[N][N];

struct coord
{
    int i, j;
};

char covered[N][N];

int coverdist[N][N];
int left_uncovered;


const int direction_count=4;
int mov[4][2]=
{
    1,0,
    0,1,
    -1,0,
    0,-1,
};

void generate_coverdist()
{
    left_uncovered=0;
    queue <coord> pq;
    for (int i=0; i<n; i++)
        for (int j=0; j<m; j++)
            if (covered[i][j])
            {
                coverdist[i][j]=0;
                pq.push(coord{i,j});
            }
            else
            {
                coverdist[i][j]=-1;
                left_uncovered++;
            }
    while (pq.size())
    {
        coord cur=pq.front();
        pq.pop();
        int i=cur.i;
        int j=cur.j;
        //printf("%d %d %d\n",i,j);
        for (int d=0; d<direction_count; d++)
        {
            int ni=i+mov[d][0];
            int nj=j+mov[d][1];
            if (coverdist[ni][nj]==-1)
            {
                coverdist[ni][nj]=coverdist[i][j]+1;
                pq.push(coord{ni,nj});
            }
        }
    }
}

vector<coord> coord_value(vector<coord> routers)
{
    bool mark[n][m];
    for (int i=0;i<n;i++)
    {
        for (int j=0;j<m;j++)
        {
            mark[i][j] = false;
        }
    }

    vector<coord> ret = {};

    for (int i=0;i<routers.size();i++)
    {
        // check if outta bounds
        if (routers[i].i < 0 || routers[i].i > n - 1 || routers[i].j < 0 || routers[i].j > m - 1) assert(false);
        // check if wallhacked
        if (board[routers[i].i][routers[i].j] == '#') assert(false);
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
                if (board[xt][yt] == '#') break; // no need to go on
                if (!mark[xt][yt])
                {
                    mark[xt][yt] = true;
                    if (board[xt][yt] == '.') ret.push_back({xt, yt});
                }
            }

            // top-right
            for (int dy=0;dy<=radius;dy++)
            {
                int yt = routers[i].j + dy;
                if (yt > m - 1) break;
                if (board[xt][yt] == '#') break;
                if (!mark[xt][yt])
                {
                    mark[xt][yt] = true;
                    if (board[xt][yt] == '.') ret.push_back({xt, yt});
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
                if (board[xt][yt] == '#') break; // no need to go on
                if (!mark[xt][yt])
                {
                    mark[xt][yt] = true;
                    if (board[xt][yt] == '.') ret.push_back({xt, yt});
                }
            }

            // top-right
            for (int dy=0;dy<=radius;dy++)
            {
                int yt = routers[i].j + dy;
                if (yt > m - 1) break;
                if (board[xt][yt] == '#') break;
                if (!mark[xt][yt])
                {
                    mark[xt][yt] = true;
                    if (board[xt][yt] == '.') ret.push_back({xt, yt});
                }
            }
        }
    }
    return ret;
}

int main()
{
    freopen(test_files[0],"r",stdin);
    //freopen("dump.txt","w",stdout);
    scanf("%d %d %d", &n, &m, &radius);
    scanf("%d %d %d", &cost_edge, &cost_router, &budget);
    scanf("%d %d", &start_i, &start_j);

    printf("%d %d\n",start_i,start_j);

    for(int i = 0; i < n; i++)
        for(int j = 0; j < m; j++)
            scanf(" %c", &board[i][j]);

    for (int i=0; i<n; i++)
        for (int j=0; j<m; j++)
        {
            if (board[i][j]!='.')
            {
                covered[i][j]=1;
            }
            else
            {
                covered[i][j]==0;
            }
        }
    vector <coord> routers;
    generate_coverdist();
    while (left_uncovered)
    {
        int ti,tj;
        for (int i=0; i<n; i++)
            for (int j=0; j<m; j++)
                if (covered[i][j]==0)
                {
                    ti=i;
                    tj=j;
                }
        vector <coord> seen=coord_value(vector<coord>{{ti,tj}});
        int soli,solj,solval;
        solval=-1;
        for (int i=0; i<seen.size(); i++)
        {
            int ni,nj;
            ni=seen[i].i;
            nj=seen[i].j;
            if (coverdist[ni][nj]>solval)
            {
                soli=ni;
                solj=nj;
                solval=coverdist[ni][nj];
            }
        }
        routers.push_back(coord{soli,solj});
        printf("%d %d\n",soli,solj);
        fflush(stdout);
        seen=coord_value(vector<coord>{{soli,solj}});
        for (int i=0; i<seen.size(); i++)
            covered[seen[i].i][seen[i].j]=1;
        generate_coverdist();
    }
    return 0;
}
