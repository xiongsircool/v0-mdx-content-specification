import { getAnnouncements, getPosts, getLearnResources, getMeetings } from "@/lib/server-markdown-loader"
import { Navigation } from "@/components/navigation"
import { Footer } from "@/components/footer"
import SponsorsSection from "@/components/sponsors-section"
import { HeroSection } from "@/components/hero-section"
import { FeaturesGrid } from "@/components/features-grid"
import { AnimatedContentSection } from "@/components/animated-content-section"
import { LiteratureShowcase } from "@/components/literature-showcase"
import { Users, BookOpen, Target, TrendingUp, ChevronRight } from "lucide-react"
import { Button } from "@/components/ui/button"
import { Badge } from "@/components/ui/badge"
import Link from "next/link"

export default async function HomePage() {
  // Get content data
  const allAnnouncements = getAnnouncements()
  const allPosts = getPosts()
  const allResources = getLearnResources()
  const allMeetings = getMeetings()

  // Calculate site statistics
  const siteStats = {
    posts: allPosts.length,
    announcements: allAnnouncements.length,
    resources: allResources.length,
    meetings: allMeetings.length
  }

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
      <HeroSection />
      <FeaturesGrid siteStats={siteStats} />
      <AnimatedContentSection announcements={recentAnnouncements} posts={recentPosts} />

      <LiteratureShowcase />

      <SponsorsSection />

      <section className="py-24 px-4 bg-gradient-to-br from-gray-900 via-blue-900 to-purple-900">
        <div className="container mx-auto text-center">
          <div className="flex flex-col items-center">
            <h2 className="text-4xl md:text-5xl font-bold mb-6 text-transparent bg-clip-text bg-gradient-to-r from-blue-400 via-purple-400 to-pink-400">
              准备好加入我们的生物信息联盟了吗？
            </h2>
            <p className="text-xl md:text-2xl text-gray-300 mb-12 max-w-3xl mx-auto leading-relaxed font-medium">
              与志同道合的伙伴一起，探索生命科学与代码的奇妙世界
            </p>
            <div className="flex flex-col sm:flex-row gap-6 justify-center">
              <div className="group">
                <Button
                  asChild
                  size="lg"
                  className="text-lg px-12 py-6 bg-gradient-to-r from-blue-600 to-purple-600 hover:from-blue-700 hover:to-purple-700 text-white shadow-2xl hover:shadow-blue-500/25 transition-all duration-300"
                >
                  <Link href="/announcements" className="flex items-center gap-3">
                    <Users className="w-6 h-6" />
                    加入联盟
                  </Link>
                </Button>
              </div>
              <div className="group">
                <Button
                  asChild
                  variant="outline"
                  size="lg"
                  className="text-lg px-12 py-6 border-2 border-blue-400 text-blue-300 hover:bg-blue-400/10 hover:text-blue-400 transition-all duration-300"
                >
                  <Link href="/posts" className="flex items-center gap-3">
                    <BookOpen className="w-6 h-6" />
                    了解分享
                  </Link>
                </Button>
              </div>
            </div>
          </div>
        </div>
      </section>

      <Footer />
    </div>
  )
}
