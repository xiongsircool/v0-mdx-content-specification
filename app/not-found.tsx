import { Navigation } from "@/components/navigation"
import { Footer } from "@/components/footer"
import { Button } from "@/components/ui/button"
import { Card, CardContent } from "@/components/ui/card"
import { ArrowLeft, Home, Search, BookOpen, Calendar, Users } from "lucide-react"
import Link from "next/link"
import { MetallicTitle } from "@/components/metallic-title"

export default function NotFound() {
  return (
    <div className="min-h-screen bg-background">
      <Navigation />

      <main className="min-h-[80vh] flex items-center justify-center px-4 py-20">
        <div className="container mx-auto max-w-4xl">
          <Card className="shadow-xl border-0 overflow-hidden">
            <CardContent className="p-0">
              <div className="relative">
                {/* Background decoration */}
                <div className="absolute inset-0 bg-gradient-to-br from-primary/5 via-accent/5 to-primary/5"></div>
                <div className="absolute top-0 right-0 w-64 h-64 bg-primary/10 rounded-full blur-3xl"></div>
                <div className="absolute bottom-0 left-0 w-48 h-48 bg-accent/10 rounded-full blur-3xl"></div>

                <div className="relative z-10 p-12 md:p-16">
                  <div className="text-center mb-8">
                    <div className="inline-flex items-center justify-center w-24 h-24 rounded-full bg-primary/10 mb-6">
                      <span className="text-6xl font-bold text-primary">404</span>
                    </div>
                    <div className="text-4xl md:text-5xl font-bold mb-4">
                      <MetallicTitle />
                    </div>
                    <p className="text-xl text-muted-foreground mb-8 max-w-2xl mx-auto">
                      抱歉，您访问的页面不存在或已被移除。
                    </p>
                  </div>

                  <div className="grid grid-cols-1 md:grid-cols-2 gap-6 mb-12">
                    <div className="space-y-4">
                      <h3 className="text-lg font-semibold flex items-center gap-2">
                        <Search className="h-5 w-5 text-primary" />
                        快速导航
                      </h3>
                      <div className="space-y-2">
                        <Button asChild variant="outline" className="w-full justify-start">
                          <Link href="/">
                            <Home className="mr-2 h-4 w-4" />
                            返回首页
                          </Link>
                        </Button>
                        <Button asChild variant="outline" className="w-full justify-start">
                          <Link href="/learn">
                            <BookOpen className="mr-2 h-4 w-4" />
                            学习资源
                          </Link>
                        </Button>
                        <Button asChild variant="outline" className="w-full justify-start">
                          <Link href="/announcements">
                            <Calendar className="mr-2 h-4 w-4" />
                            最新公告
                          </Link>
                        </Button>
                        <Button asChild variant="outline" className="w-full justify-start">
                          <Link href="/posts">
                            <BookOpen className="mr-2 h-4 w-4" />
                            技术文章
                          </Link>
                        </Button>
                      </div>
                    </div>

                    <div className="space-y-4">
                      <h3 className="text-lg font-semibold flex items-center gap-2">
                        <Users className="h-5 w-5 text-primary" />
                        常用链接
                      </h3>
                      <div className="space-y-2">
                        <Button asChild variant="outline" className="w-full justify-start">
                          <Link href="/meetings">
                            <Calendar className="mr-2 h-4 w-4" />
                            会议活动
                          </Link>
                        </Button>
                        <Button asChild variant="outline" className="w-full justify-start">
                          <Link href="/about">
                            <Users className="mr-2 h-4 w-4" />
                            关于我们
                          </Link>
                        </Button>
                        <Button asChild variant="outline" className="w-full justify-start">
                          <Link href="/search">
                            <Search className="mr-2 h-4 w-4" />
                            全站搜索
                          </Link>
                        </Button>
                      </div>
                    </div>
                  </div>

                  <div className="text-center">
                    <p className="text-muted-foreground mb-6">
                      如果您认为这是一个错误，请尝试：
                    </p>
                    <div className="flex flex-wrap justify-center gap-4 text-sm">
                      <span className="px-3 py-1 bg-primary/10 rounded-full text-primary">
                        检查URL拼写
                      </span>
                      <span className="px-3 py-1 bg-primary/10 rounded-full text-primary">
                        使用搜索功能
                      </span>
                      <span className="px-3 py-1 bg-primary/10 rounded-full text-primary">
                        从首页重新导航
                      </span>
                    </div>
                  </div>

                  <div className="text-center mt-12">
                    <Button asChild size="lg" className="premium-gradient">
                      <Link href="/">
                        <ArrowLeft className="mr-2 h-5 w-5" />
                        返回首页
                      </Link>
                    </Button>
                  </div>
                </div>
              </div>
            </CardContent>
          </Card>
        </div>
      </main>

      <Footer />
    </div>
  )
}