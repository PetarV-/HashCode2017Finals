#include <vector>
#include <string>
#include <queue>
#include <cstdio>
#include <cassert>
#include <algorithm>

using namespace std;

char test_files[4][50]={"tests/charleston_road.in","tests/lets_go_higher.in","tests/opera.in","tests/rue_de_londres.in"};

const int N = 1000;
int n, m, radius, cost_edge, cost_router, budget, start_i, start_j;
char board[N][N];

struct coord
{
    int i, j;
    inline bool operator != (const coord &op) const
    {
        return i!=op.i || j!=op.j;
    }
    inline bool operator < (const coord &op) const
    {
        if (i!=op.i)
            return i<op.i;
        return j<op.j;
    }
};

char covered[N][N];

int coverdist[N][N];
int left_uncovered;


const int direction_count=8;
int mov[8][2]=
{
    1,0,
    0,1,
    -1,0,
    0,-1,
    1,1,
    1,-1,
    -1,1,
    -1,-1,
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

        int ly = -1;
        int ry = m;

        for (int dx=0;dx>=-radius;dx--)
        {
            int xt = routers[i].i + dx;
            if (xt < 0) break;
            if (board[xt][routers[i].j] == '#') break;

            // top-left
            for (int dy=0;dy>=-radius;dy--)
            {
                int yt = routers[i].j + dy;
                if (yt <= ly) break;
                if (board[xt][yt] == '#')
                {
                    ly = max(ly, yt);
                    break; // no need to go on
                }
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
                if (yt >= ry) break;
                if (board[xt][yt] == '#')
                {
                    ry = min(ry, yt);
                    break;
                }
                if (!mark[xt][yt])
                {
                    mark[xt][yt] = true;
                    if (board[xt][yt] == '.') ret.push_back({xt, yt});
                }
            }
        }

        ly = -1;
        ry = m;

        // now expand down
        for (int dx=0;dx<=radius;dx++)
        {
            int xt = routers[i].i + dx;
            if (xt > n - 1) break;
            if (board[xt][routers[i].j] == '#') break;

            // bottom-left
            for (int dy=0;dy>=-radius;dy--)
            {
                int yt = routers[i].j + dy;
                if (yt <= ly) break;
                if (board[xt][yt] == '#')
                {
                    ly = max(ly, yt);
                    break; // no need to go on
                }
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
                if (yt >= ry) break;
                if (board[xt][yt] == '#')
                {
                    ry = min(ry, yt);
                    break;
                }
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

// -----------
// output

void write_output(string filename, vector<coord> routers, vector<coord> backbone)
{
    FILE *f = fopen(filename.c_str(), "w");
    if(!f)
    {
		fprintf(stderr, "Opening file %s for writing failed!\n", filename.c_str());
		return;
    }

    fprintf(f, "%d\n", backbone.size());
	for(coord i : backbone)
		fprintf(f, "%d %d\n", i.i, i.j);

	fprintf(f, "%d\n", routers.size());
	for(coord i : routers)
		fprintf(f, "%d %d\n", i.i, i.j);

	fclose(f);
}

// ----------
// ----------
// backbone

struct edge
{
	int from, to;
	int len;
} ;

bool operator<(edge a, edge b)
{
	// > for priority_queue as min-heap
	return a.len == b.len ? a.from == b.from ? a.to < b.to : a.from < b.from : a.len > b.len;
}

int dist(coord a, coord b)
{
	return max(abs(a.i - b.i), abs(a.j - b.j));
}

int sgn(int x) { return x ? (x > 0 ? 1 : -1) : 0; }

bool in_mst[N];
bool backbone[N][N];
int backcomp[N][N];
coord backedge[N][N];

vector <coord> markcomp(coord start,int col)
{
    queue <coord> pq;
    vector <coord> comp_parts;
    pq.push(start);
    backcomp[start.i][start.j]=col;
    backedge[start.i][start.j]=coord{start.i,start.j};
    while (pq.size())
    {
        coord cur=pq.front();
        pq.pop();
        comp_parts.push_back(cur);
        for (int d=0; d<direction_count; d++)
        {
            int ni=cur.i+mov[d][0];
            int nj=cur.j+mov[d][1];
            if (backcomp[ni][nj]==0 && backbone[ni][nj]==true)
            {
                pq.push(coord{ni,nj});
                backcomp[ni][nj]=col;
                backedge[ni][nj]=coord{ni,nj};
            }
        }
    }
    return comp_parts;
}

vector<coord> make_backbone(vector<coord> routers)
{
    backbone[start_i][start_j]=true;
    for (int i=0; i<routers.size(); i++)
        backbone[routers[i].i][routers[i].j]=true;
    int compnum;
    do
    {
        compnum=0;
        for (int i=0; i<n; i++)
            for (int j=0; j<m; j++)
            {
                backcomp[i][j]=0;
                backedge[i][j]=coord{0,0};
            }
        queue <coord> pq;
        for (int i=0; i<n; i++)
            for (int j=0; j<m; j++)
                if (backbone[i][j] && backcomp[i][j]==0)
                {
                    compnum++;
                    vector <coord> comp=markcomp(coord{i,j},compnum);
                    for (int i=0; i<comp.size(); i++)
                        pq.push(comp[i]);
                }
        printf("%d\n",compnum);
        if (compnum>1)
        {
            while (pq.size())
            {
                coord cur=pq.front();
                //printf("PQ SIZE %d\n",pq.size());
                //printf("POINT %d %d\n",cur.i,cur.j);
                pq.pop();
                for (int d=0; d<direction_count; d++)
                {
                    int ni=cur.i+mov[d][0];
                    int nj=cur.j+mov[d][1];
                    if (backcomp[ni][nj]==0)
                    {
                        backcomp[ni][nj]=backcomp[cur.i][cur.j];
                        backedge[ni][nj]=cur;
                        pq.push(coord{ni,nj});
                        continue;
                    }
                    if (backcomp[ni][nj]!=backcomp[cur.i][cur.j])
                    {
                        printf("COMPS %d %d\n",backcomp[ni][nj],backcomp[cur.i][cur.j]);
                        coord addpos=coord{ni,nj};
                        while (backedge[addpos.i][addpos.j]!=addpos)
                        {
                            backbone[addpos.i][addpos.j]=true;
                            addpos=backedge[addpos.i][addpos.j];
                        }
                        addpos=cur;
                        while (backedge[addpos.i][addpos.j]!=addpos)
                        {
                            backbone[addpos.i][addpos.j]=true;
                            addpos=backedge[addpos.i][addpos.j];
                        }
                        while (pq.size())
                            pq.pop();
                        break;
                    }
                }
            }
        }
    } while (compnum>1);
    vector <coord> backbone_fields;
    backbone_fields.push_back(coord{start_i,start_j});
    backbone[start_i][start_j]=0;
    for (int i=0; i<backbone_fields.size(); i++)
        for (int d=0; d<direction_count; d++)
        {
            int ni=backbone_fields[i].i+mov[d][0];
            int nj=backbone_fields[i].j+mov[d][1];
            if (backbone[ni][nj])
            {
                backbone[ni][nj]=0;
                backbone_fields.push_back(coord{ni,nj});
            }
        }
    backbone_fields.erase(backbone_fields.begin());
    return backbone_fields;
}

// ----------


int main()
{
    freopen(test_files[1],"r",stdin);
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
    int lastuncoveredi=0;
    int bindwall=0;
    while (left_uncovered)
    {
        if (routers.size()%50==0)
            printf("PROGRESS %d %d\n",left_uncovered,routers.size());
        int ti,tj,bestwall;
        bestwall=-1;
        for (int i=0; i<n; i++)
        {
            for (int j=0; j<m; j++)
                if (covered[i][j]==0)
                    if (coverdist[i][j]>bestwall)
                    {
                        ti=i;
                        tj=j;
                        bestwall=coverdist[i][j];
                    }
        }
        int soli,solj;
        if (bestwall>3)
        {
            soli=ti;
            solj=tj;
        }
        else
        {
            bindwall=1;
            break;
        }
        routers.push_back(coord{soli,solj});
        vector <coord> seen=coord_value(vector<coord>{{soli,solj}});
        for (int i=0; i<seen.size(); i++)
            covered[seen[i].i][seen[i].j]=1;
        generate_coverdist();
    }
    priority_queue <pair <int,coord> > fields_covered;
    for (int i=0; i<n; i++)
        for (int j=0; j<m; j++)
            if (board[i][j]=='.')
            {
                vector <coord> sees=coord_value(vector<coord>{{i,j}});
                int cscore=0;
                for (int t=0; t<sees.size(); t++)
                    if (covered[sees[t].i][sees[t].j]==0)
                        cscore++;
                fields_covered.push(make_pair(cscore,coord{i,j}));
            }
    while (left_uncovered)
    {
         pair <int,coord> najbolji=fields_covered.top();
        fields_covered.pop();
        int i=najbolji.second.i;
        int j=najbolji.second.j;

        vector <coord> sees=coord_value(vector<coord>{{i,j}});
        int cscore=0;
        for (int t=0; t<sees.size(); t++)
            if (covered[sees[t].i][sees[t].j]==0)
                cscore++;
        int soli,solj;
        bool passed;
        if (cscore==najbolji.first)
        {
            soli=i;
            solj=j;
            passed=true;
        }
        else
        {
            fields_covered.push(make_pair(cscore,najbolji.second));
            passed=false;
        }

        if (passed)
        {
            routers.push_back(coord{soli,solj});
            vector <coord> seen=coord_value(vector<coord>{{soli,solj}});
            for (int i=0; i<seen.size(); i++)
            {
                if (covered[seen[i].i][seen[i].j]==0)
                {
                    covered[seen[i].i][seen[i].j]=1;
                    left_uncovered--;
                }
            }
        }
    }

    write_output("out/charleston_road.out", routers, make_backbone(routers));
    return 0;
}
