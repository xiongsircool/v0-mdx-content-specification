import { allAnnouncements, allPosts } from "@/lib/contentlayer-mock"
import { AnnouncementCard, PostCard } from "@/components/content-card"
import { Navigation } from "@/components/navigation"
import { Footer } from "@/components/footer"
import { Button } from "@/components/ui/button"
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card"
import { Badge } from "@/components/ui/badge"
import { Megaphone, BookOpen, Users, TrendingUp, ArrowRight, Calendar, Target } from "lucide-react"
import Link from "next/link"
import { MetallicTitle } from "@/components/metallic-title"

export default function HomePage() {
  // Get recent announcements (non-expired, sorted by priority and date)
  const recentAnnouncements = allAnnouncements
    .filter((announcement) => !announcement.isExpired)
    .sort((a, b) => {
      if (a.priority !== b.priority) {
        return b.priority - a.priority
      }
      return new Date(b.publishedAt).getTime() - new Date(a.publishedAt).getTime()
    })
    .slice(0, 3)

  // Get recent posts
  const recentPosts = allPosts
    .sort((a, b) => new Date(b.publishedAt).getTime() - new Date(a.publishedAt).getTime())
    .slice(0, 3)

  return (
    <div className="min-h-screen bg-background">
      <Navigation />

      <section className="relative py-20 px-4 animated-bg min-h-[80vh] flex items-center">
        {/* Floating elements for animation */}
        <div className="floating-element"></div>
        <div className="floating-element"></div>
        <div className="floating-element"></div>
        <div className="floating-element"></div>
        <div className="glow-orb"></div>
        <div className="glow-orb"></div>

        <div className="container mx-auto text-center relative z-10">
          <div className="flex justify-center mb-8">
            <div className="flex h-20 w-20 items-center justify-center rounded-3xl premium-gradient shadow-2xl">
              <Target className="h-10 w-10 text-primary" />
            </div>
          </div>
          <MetallicTitle />
          <p className="text-xl md:text-2xl text-gray-700 max-w-4xl mx-auto mb-12 text-pretty leading-relaxed font-medium">
            探索生命科学与计算技术的交汇点，培养下一代生物信息学人才。
            我们致力于提供系统化的学习资源、实践机会和学术交流平台。
          </p>
          <div className="flex flex-col sm:flex-row gap-6 justify-center">
            <Button
              asChild
              size="lg"
              className="text-lg px-8 py-6 shadow-lg hover:shadow-xl transition-all duration-300 premium-gradient"
            >
              <Link href="/learn">
                开始学习
                <ArrowRight className="ml-2 h-5 w-5" />
              </Link>
            </Button>
            <Button
              asChild
              variant="outline"
              size="lg"
              className="text-lg px-8 py-6 premium-gradient shadow-lg hover:shadow-xl transition-all duration-300 bg-transparent"
            >
              <Link href="/about">了解更多</Link>
            </Button>
          </div>
        </div>
      </section>

      <section className="py-20 px-4 bg-gradient-to-b from-white to-secondary/30">
        <div className="container mx-auto">
          <div className="grid grid-cols-1 md:grid-cols-3 gap-8">
            <Card className="shadow-lg hover:shadow-xl transition-all duration-300 border-0 premium-gradient">
              <CardContent className="p-8 text-center">
                <div className="flex justify-center mb-6">
                  <div className="flex h-16 w-16 items-center justify-center rounded-2xl bg-gradient-to-br from-blue-500/20 to-blue-600/20 border border-blue-500/20">
                    <Users className="h-8 w-8 text-blue-600" />
                  </div>
                </div>
                <h3 className="text-3xl font-bold mb-3 text-foreground">50+</h3>
                <p className="text-muted-foreground text-lg">活跃成员</p>
              </CardContent>
            </Card>
            <Card className="shadow-lg hover:shadow-xl transition-all duration-300 border-0 premium-gradient">
              <CardContent className="p-8 text-center">
                <div className="flex justify-center mb-6">
                  <div className="flex h-16 w-16 items-center justify-center rounded-2xl bg-gradient-to-br from-green-500/20 to-green-600/20 border border-green-500/20">
                    <BookOpen className="h-8 w-8 text-green-600" />
                  </div>
                </div>
                <h3 className="text-3xl font-bold mb-3 text-foreground">{allPosts.length}</h3>
                <p className="text-muted-foreground text-lg">技术文章</p>
              </CardContent>
            </Card>
            <Card className="shadow-lg hover:shadow-xl transition-all duration-300 border-0 premium-gradient">
              <CardContent className="p-8 text-center">
                <div className="flex justify-center mb-6">
                  <div className="flex h-16 w-16 items-center justify-center rounded-2xl bg-gradient-to-br from-purple-500/20 to-purple-600/20 border border-purple-500/20">
                    <TrendingUp className="h-8 w-8 text-purple-600" />
                  </div>
                </div>
                <h3 className="text-3xl font-bold mb-3 text-foreground">10+</h3>
                <p className="text-muted-foreground text-lg">实践项目</p>
              </CardContent>
            </Card>
          </div>
        </div>
      </section>

      {/* Recent Announcements */}
      <section className="py-16 px-4">
        <div className="container mx-auto">
          <div className="flex items-center justify-between mb-8">
            <div className="flex items-center gap-3">
              <div className="flex h-10 w-10 items-center justify-center rounded-lg bg-primary/10">
                <Megaphone className="h-5 w-5 text-primary" />
              </div>
              <div>
                <h2 className="text-3xl font-bold">最新公告</h2>
                <p className="text-muted-foreground">了解俱乐部最新动态</p>
              </div>
            </div>
            <Button asChild variant="outline">
              <Link href="/announcements">查看全部</Link>
            </Button>
          </div>

          {recentAnnouncements.length > 0 ? (
            <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
              {recentAnnouncements.map((announcement) => (
                <AnnouncementCard key={announcement.slug} announcement={announcement} />
              ))}
            </div>
          ) : (
            <Card>
              <CardContent className="p-12 text-center">
                <Megaphone className="h-16 w-16 text-muted-foreground mx-auto mb-4" />
                <h3 className="text-xl font-semibold mb-2">暂无公告</h3>
                <p className="text-muted-foreground">目前没有发布任何公告</p>
              </CardContent>
            </Card>
          )}
        </div>
      </section>

      {/* Recent Posts */}
      <section className="py-16 px-4 bg-muted/50">
        <div className="container mx-auto">
          <div className="flex items-center justify-between mb-8">
            <div className="flex items-center gap-3">
              <div className="flex h-10 w-10 items-center justify-center rounded-lg bg-primary/10">
                <BookOpen className="h-5 w-5 text-primary" />
              </div>
              <div>
                <h2 className="text-3xl font-bold">最新文章</h2>
                <p className="text-muted-foreground">深度技术分享与教程</p>
              </div>
            </div>
            <Button asChild variant="outline">
              <Link href="/posts">查看全部</Link>
            </Button>
          </div>

          {recentPosts.length > 0 ? (
            <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
              {recentPosts.map((post) => (
                <PostCard key={post.slug} post={post} />
              ))}
            </div>
          ) : (
            <Card>
              <CardContent className="p-12 text-center">
                <BookOpen className="h-16 w-16 text-muted-foreground mx-auto mb-4" />
                <h3 className="text-xl font-semibold mb-2">暂无文章</h3>
                <p className="text-muted-foreground">目前没有发布任何技术文章</p>
              </CardContent>
            </Card>
          )}
        </div>
      </section>

      {/* Learning Tracks */}
      <section className="py-16 px-4">
        <div className="container mx-auto">
          <div className="text-center mb-12">
            <h2 className="text-3xl font-bold mb-4">学习路径</h2>
            <p className="text-xl text-gray-600 max-w-2xl mx-auto font-medium">
              系统化的学习资源，帮助你从零基础到专业水平
            </p>
          </div>

          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
            <Card className="transition-all hover:shadow-md">
              <CardHeader>
                <div className="flex items-center gap-3">
                  <div className="flex h-10 w-10 items-center justify-center rounded-lg bg-blue-100 dark:bg-blue-900">
                    <Target className="h-5 w-5 text-blue-600 dark:text-blue-400" />
                  </div>
                  <CardTitle>生物信息学基础</CardTitle>
                </div>
              </CardHeader>
              <CardContent>
                <p className="text-muted-foreground mb-4">从生物学概念到计算方法的全面介绍</p>
                <div className="flex items-center gap-2 mb-4">
                  <Badge variant="outline">入门友好</Badge>
                  <Badge variant="secondary">理论基础</Badge>
                </div>
                <Button asChild variant="outline" size="sm" className="w-full bg-transparent">
                  <Link href="/learn">开始学习</Link>
                </Button>
              </CardContent>
            </Card>

            <Card className="transition-all hover:shadow-md">
              <CardHeader>
                <div className="flex items-center gap-3">
                  <div className="flex h-10 w-10 items-center justify-center rounded-lg bg-green-100 dark:bg-green-900">
                    <BookOpen className="h-5 w-5 text-green-600 dark:text-green-400" />
                  </div>
                  <CardTitle>编程技能</CardTitle>
                </div>
              </CardHeader>
              <CardContent>
                <p className="text-muted-foreground mb-4">Python、R语言等生物信息学编程技能</p>
                <div className="flex items-center gap-2 mb-4">
                  <Badge variant="outline">实践导向</Badge>
                  <Badge variant="secondary">项目驱动</Badge>
                </div>
                <Button asChild variant="outline" size="sm" className="w-full bg-transparent">
                  <Link href="/learn">开始学习</Link>
                </Button>
              </CardContent>
            </Card>

            <Card className="transition-all hover:shadow-md">
              <CardHeader>
                <div className="flex items-center gap-3">
                  <div className="flex h-10 w-10 items-center justify-center rounded-lg bg-purple-100 dark:bg-purple-900">
                    <TrendingUp className="h-5 w-5 text-purple-600 dark:text-purple-400" />
                  </div>
                  <CardTitle>单细胞分析</CardTitle>
                </div>
              </CardHeader>
              <CardContent>
                <p className="text-muted-foreground mb-4">单细胞RNA测序数据分析的前沿技术</p>
                <div className="flex items-center gap-2 mb-4">
                  <Badge variant="outline">前沿技术</Badge>
                  <Badge variant="secondary">高级应用</Badge>
                </div>
                <Button asChild variant="outline" size="sm" className="w-full bg-transparent">
                  <Link href="/learn">开始学习</Link>
                </Button>
              </CardContent>
            </Card>
          </div>
        </div>
      </section>

      <section className="py-24 px-4 bg-gradient-to-r from-primary/5 via-accent/5 to-primary/5">
        <div className="container mx-auto text-center">
          <h2 className="text-4xl md:text-5xl font-bold mb-6 text-balance">准备好开始你的生物信息学之旅了吗？</h2>
          <p className="text-xl md:text-2xl text-gray-600 mb-12 max-w-3xl mx-auto text-pretty leading-relaxed font-medium">
            加入我们的社区，获得专业指导、实践机会和同行交流
          </p>
          <div className="flex flex-col sm:flex-row gap-6 justify-center">
            <Button
              asChild
              size="lg"
              className="text-lg px-8 py-6 shadow-lg hover:shadow-xl transition-all duration-300 premium-gradient"
            >
              <Link href="/learn">
                <Calendar className="mr-2 h-5 w-5" />
                开始学习
              </Link>
            </Button>
            <Button
              asChild
              variant="outline"
              size="lg"
              className="text-lg px-8 py-6 premium-gradient shadow-lg hover:shadow-xl transition-all duration-300 bg-transparent"
            >
              <Link href="/announcements">查看招新信息</Link>
            </Button>
          </div>
        </div>
      </section>

      <Footer />
    </div>
  )
}
